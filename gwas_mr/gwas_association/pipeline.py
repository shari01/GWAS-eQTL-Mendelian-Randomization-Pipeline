"""End-to-end GWAS to transcriptomics pipeline orchestrator."""

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Set, Union

import pandas as pd
from openai import OpenAI

from .biosample import (
    enforce_disease_tissue_consistency,
    llm_filter_tissue,
    normalize_biosample_input,
)
from .config import MAX_STUDY_ROWS_FOR_LLM, MODEL_NAME, get_openai_client
from .filtering import (
    choose_best_preview_traits,
    fallback_rank_candidates,
    filter_rows_by_preview_traits,
    get_all_disease_values,
    preview_closest_traits,
    simple_trait_filter,
    strict_disease_filter,
)
from .genes import (
    build_gene_set_from_input,
    expand_mapped_rows_per_gene,
    find_gene_column,
    map_input_genes_to_gwas,
)
from .helpers import clean_name, get_trait_match_column, normalize_text
from .io_utils import read_table_auto
from .llm_utils import (
    correct_disease_name,
    generate_disease_variants,
    llm_rank_rows,
)
from .rag import build_or_load_gwas_rag_index
from .reporting import make_knowledge_graph


@dataclass
class PipelineResult:
    """Outputs produced by a single GWAS pipeline run."""

    run_dir: Path
    genotype_association_file: Path
    patient_transcriptomics_file: Path
    patient_graph_file: Path
    genotype_graph_file: Path
    disease_used: str
    biosample_used: str
    gwas_rows_selected: int
    genes_mapped: int


def _init_logger() -> logging.Logger:
    """Return the ``gwas_association`` logger, bootstrapping it on first call."""
    logger = logging.getLogger("gwas_association")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-7s | %(message)s",
            datefmt="%H:%M:%S",
        ))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    return logger


def run_gwas_association(
    gwas_index_file: Path,
    input_genes: Union[List[str], Sequence[str], Path],
    output_dir: Path,
    disease: str,
    biosample: str,
    client: Optional[OpenAI] = None,
    use_simple_filter: bool = True,
    model_name: str = MODEL_NAME,
    max_study_rows_for_llm: int = MAX_STUDY_ROWS_FOR_LLM,
    rag_cache_dir: Optional[Path] = None,
    max_workers: int = 10,
) -> PipelineResult:
    """Run the full GWAS to transcriptomics mapping pipeline.

    Args:
        gwas_index_file: Path to the GWAS association TSV.
        input_genes: Gene input — either a **list of gene symbols**
            (e.g. ``["FLT3", "NPM1", "RUNX1"]``) or a **Path** to a
            tabular file whose first (or auto-detected) column contains
            gene symbols.
        output_dir: Root directory for all pipeline outputs.
        disease: Disease name (will be spell-corrected automatically).
        biosample: Tissue / biosample term (will be normalised).
        client: OpenAI client.  Created via :func:`get_openai_client`
            when *None*.
        use_simple_filter: When *True*, use the faster keyword + LLM tissue
            filter.  When *False*, use the full strict + LLM ranking path.
        model_name: OpenAI model identifier.
        max_study_rows_for_llm: Cap on candidate rows sent to LLM ranking
            (strict-filter mode only).
        rag_cache_dir: Override for the RAG index cache directory.
        max_workers: Number of concurrent threads for LLM API calls.
            Higher values speed up tissue filtering and row ranking but
            may hit API rate limits.

    Returns:
        A :class:`PipelineResult` with output paths and summary statistics.

    Raises:
        FileNotFoundError: If required input files are missing.
        ValueError: If the GWAS index is missing required columns or no
            candidates match the disease query.
    """
    logger = _init_logger()

    if client is None:
        client = get_openai_client()

    t_start = time.perf_counter()

    # ==================================================================
    # Step 1: Load GWAS index
    # ==================================================================
    if not gwas_index_file.exists():
        raise FileNotFoundError(
            f"GWAS index file not found: {gwas_index_file}"
        )

    logger.info("[1/9] Loading GWAS index file...")
    gwas_df = pd.read_csv(gwas_index_file, sep="\t", low_memory=False)
    logger.info("Loaded GWAS rows: %s", f"{len(gwas_df):,}")

    required_cols = ["DISEASE/TRAIT", "STUDY", "REPORTED GENE(S)"]
    missing = [c for c in required_cols if c not in gwas_df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in GWAS index: {missing}"
        )

    # ==================================================================
    # Step 2: Build / load RAG index
    # ==================================================================
    logger.info("[2/9] Building / loading RAG index...")
    gwas_rag_index = build_or_load_gwas_rag_index(
        gwas_df,
        gwas_index_file,
        rag_cache_dir=rag_cache_dir,
        logger=logger,
    )
    trait_match_col = gwas_rag_index.get(
        "trait_column", get_trait_match_column(gwas_df)
    )
    logger.info("Trait column for matching: %s", trait_match_col)

    # ==================================================================
    # Step 3: Resolve input genes
    # ==================================================================
    if isinstance(input_genes, Path):
        if not input_genes.exists():
            raise FileNotFoundError(f"Input file not found: {input_genes}")
        logger.info("[3/9] Loading gene file: %s", input_genes.name)
        input_df = read_table_auto(input_genes)
        try:
            gene_col = find_gene_column(input_df)
        except ValueError:
            gene_col = input_df.columns[0]
            logger.info("No standard gene column found; using first column: %s", gene_col)
        input_gene_set: Set[str] = build_gene_set_from_input(input_df, gene_col)
        logger.info("Detected gene column: %s", gene_col)
    else:
        logger.info("[3/9] Received gene list directly (%d entries)", len(input_genes))
        input_gene_set = {g.strip().upper() for g in input_genes if g.strip()}

    if not input_gene_set:
        raise ValueError("No valid genes provided in input_genes.")
    logger.info("Unique input genes: %s", f"{len(input_gene_set):,}")

    # ==================================================================
    # Step 4: Normalise biosample & correct disease spelling
    # ==================================================================
    logger.info("[4/9] Normalising disease and biosample terms...")
    normalized_biosample = normalize_biosample_input(
        biosample,
        client=client,
        model_name=model_name,
        logger=logger,
    )
    if normalize_text(normalized_biosample) != normalize_text(biosample):
        logger.info(
            "Biosample normalised: '%s' -> '%s'",
            biosample, normalized_biosample,
        )

    all_disease_values = get_all_disease_values(gwas_df)
    corrected_disease = correct_disease_name(
        disease,
        all_disease_values,
        client=client,
        model_name=model_name,
        logger=logger,
    )
    disease_variants = [corrected_disease]

    if normalize_text(corrected_disease) != normalize_text(disease):
        logger.info(
            "Disease corrected: '%s' -> '%s'", disease, corrected_disease,
        )

    # ==================================================================
    # Step 5 + 6: Disease filtering & tissue selection
    # ==================================================================
    if use_simple_filter:
        llm_selected_df, ranked_df = _run_simple_filter(
            gwas_df=gwas_df,
            corrected_disease=corrected_disease,
            normalized_biosample=normalized_biosample,
            trait_match_col=trait_match_col,
            client=client,
            model_name=model_name,
            logger=logger,
            max_workers=max_workers,
        )
    else:
        llm_selected_df, ranked_df, disease_variants = _run_strict_filter(
            gwas_df=gwas_df,
            corrected_disease=corrected_disease,
            normalized_biosample=normalized_biosample,
            trait_match_col=trait_match_col,
            gwas_rag_index=gwas_rag_index,
            client=client,
            model_name=model_name,
            max_study_rows_for_llm=max_study_rows_for_llm,
            logger=logger,
            max_workers=max_workers,
        )

    # ==================================================================
    # Step 7: Map patient genes to GWAS rows
    # ==================================================================
    logger.info("[7/9] Mapping input genes to GWAS rows...")
    mapped_df = map_input_genes_to_gwas(llm_selected_df, input_gene_set)
    mapped_expanded_df = expand_mapped_rows_per_gene(mapped_df)
    logger.info("Mapped GWAS rows: %s", f"{len(mapped_df):,}")
    logger.info("Expanded per-gene rows: %s", f"{len(mapped_expanded_df):,}")

    # ==================================================================
    # Step 8: Save output files
    # ==================================================================
    logger.info("[8/9] Saving outputs...")
    disease_tag = clean_name(corrected_disease)
    biosample_tag = clean_name(normalized_biosample)
    run_name = f"{disease_tag}_{biosample_tag}"
    run_dir = output_dir / "gwas_association"
    run_dir.mkdir(parents=True, exist_ok=True)

    population_out = run_dir / f"{run_name}_gwas_genotype_association.csv"
    transcriptomics_out = run_dir / f"{run_name}_gwas_patient_transcriptomics.csv"
    patient_graph_out = run_dir / f"{run_name}_patient_graph.html"
    genotype_graph_out = run_dir / f"{run_name}_genotype_graph.html"

    llm_selected_df.to_csv(population_out, index=False)
    mapped_df.to_csv(transcriptomics_out, index=False)
    logger.info("Saved: %s", population_out)
    logger.info("Saved: %s", transcriptomics_out)

    # ==================================================================
    # Step 9: Generate knowledge graphs
    # ==================================================================
    logger.info("[9/9] Generating knowledge graphs...")
    try:
        make_knowledge_graph(
            disease_name=corrected_disease,
            selected_df=mapped_df,
            mapped_expanded_df=mapped_expanded_df,
            output_path=patient_graph_out,
            logger=logger,
        )
        make_knowledge_graph(
            disease_name=corrected_disease,
            selected_df=llm_selected_df,
            mapped_expanded_df=mapped_expanded_df,
            output_path=genotype_graph_out,
            logger=logger,
        )
    except Exception as e:
        logger.warning("Graph generation failed: %s", e)
        logger.warning("Continuing pipeline without graphs.")

    # ==================================================================
    # Done
    # ==================================================================
    elapsed = time.perf_counter() - t_start
    minutes, seconds = divmod(elapsed, 60)
    logger.info("Pipeline finished.  Run directory: %s", run_dir)
    logger.info(
        "Total execution time: %dm %.1fs (%.1fs)",
        int(minutes), seconds, elapsed,
    )

    return PipelineResult(
        run_dir=run_dir,
        genotype_association_file=population_out,
        patient_transcriptomics_file=transcriptomics_out,
        patient_graph_file=patient_graph_out,
        genotype_graph_file=genotype_graph_out,
        disease_used=corrected_disease,
        biosample_used=normalized_biosample,
        gwas_rows_selected=len(llm_selected_df),
        genes_mapped=len(mapped_df),
    )


# ── Internal helpers (keep the main function readable) ───────────────────


def _run_simple_filter(
    gwas_df: pd.DataFrame,
    corrected_disease: str,
    normalized_biosample: str,
    trait_match_col: str,
    client: OpenAI,
    model_name: str,
    logger: logging.Logger,
    max_workers: int = 10,
) -> tuple:
    """Simple-mode path: keyword match + LLM tissue check.

    Returns:
        ``(llm_selected_df, ranked_df)`` tuple.
    """
    logger.info("[5/9] Simple filter: keyword match on %s...", trait_match_col)

    candidate_df = simple_trait_filter(
        gwas_df, corrected_disease, trait_match_col,
    )
    if candidate_df.empty:
        raise ValueError(
            f"No GWAS rows found where {trait_match_col} "
            f"matches '{corrected_disease}'."
        )

    candidate_df = candidate_df.drop_duplicates().reset_index(drop=True)
    logger.info("Keyword-matched rows: %s", f"{len(candidate_df):,}")

    logger.info("[6/9] LLM tissue selection...")
    candidate_df = enforce_disease_tissue_consistency(
        candidate_df, corrected_disease, normalized_biosample, logger=logger,
    )
    if candidate_df.empty:
        raise ValueError(
            f"Disease-tissue mismatch: "
            f"'{corrected_disease}' vs '{normalized_biosample}'."
        )

    tissue_filtered_df = llm_filter_tissue(
        candidate_df,
        corrected_disease,
        normalized_biosample,
        client=client,
        model_name=model_name,
        logger=logger,
        max_workers=max_workers,
    )

    if tissue_filtered_df.empty:
        logger.warning(
            "LLM found no tissue match — keeping disease-only rows."
        )
        llm_selected_df = candidate_df.copy()
    else:
        llm_selected_df = tissue_filtered_df.copy()

    logger.info(
        "Rows after tissue selection: %s", f"{len(llm_selected_df):,}"
    )
    return llm_selected_df, llm_selected_df.copy()


def _run_strict_filter(
    gwas_df: pd.DataFrame,
    corrected_disease: str,
    normalized_biosample: str,
    trait_match_col: str,
    gwas_rag_index: dict,
    client: OpenAI,
    model_name: str,
    max_study_rows_for_llm: int,
    logger: logging.Logger,
    max_workers: int = 10,
) -> tuple:
    """Strict-mode path: full filtering + LLM ranking.

    Returns:
        ``(llm_selected_df, ranked_df, disease_variants)`` tuple.
    """
    logger.info("[5/9] Strict disease filtering...")
    disease_variants = [corrected_disease]

    preview_df = preview_closest_traits(
        gwas_df, corrected_disease, rag_index=gwas_rag_index, top_n=15,
    )

    if not preview_df.empty:
        for _, row in preview_df.head(5).iterrows():
            logger.info(
                "  score=%.3f | rows=%d | trait=%s",
                row["SIMILARITY_SCORE"],
                int(row["ROW_COUNT"]),
                row["DISEASE_TRAIT"],
            )

    selected_preview_traits = choose_best_preview_traits(
        preview_df, min_score=0.80, max_traits=8,
    )

    candidate_df = strict_disease_filter(
        gwas_df, corrected_disease, trait_match_col,
    )
    candidate_df = filter_rows_by_preview_traits(
        candidate_df, selected_preview_traits, trait_match_col,
    )

    if candidate_df.empty:
        logger.info("No strict match — generating disease variants...")
        disease_variants = generate_disease_variants(
            corrected_disease,
            client=client,
            model_name=model_name,
            logger=logger,
        )
        variant_frames = [
            strict_disease_filter(gwas_df, v, trait_match_col)
            for v in disease_variants
        ]
        variant_frames = [df for df in variant_frames if not df.empty]
        if variant_frames:
            candidate_df = pd.concat(variant_frames, ignore_index=True)
            candidate_df = filter_rows_by_preview_traits(
                candidate_df, selected_preview_traits, trait_match_col,
            )

    if candidate_df.empty:
        logger.info("Using fallback similarity ranking...")
        candidate_df = fallback_rank_candidates(
            gwas_df, corrected_disease, top_n=200,
        )

    if candidate_df.empty:
        raise ValueError(
            "No candidate GWAS rows found for the disease query."
        )

    candidate_df = candidate_df.drop_duplicates().reset_index(drop=True)
    logger.info("Disease candidate rows: %s", f"{len(candidate_df):,}")

    if len(candidate_df) > max_study_rows_for_llm:
        logger.info(
            "Capping candidates to %d before LLM ranking.",
            max_study_rows_for_llm,
        )
        candidate_df = candidate_df.head(max_study_rows_for_llm).copy()

    logger.info("[6/9] LLM ranking by disease + tissue + study...")
    ranked_df = llm_rank_rows(
        candidate_df,
        corrected_disease,
        normalized_biosample,
        client=client,
        model_name=model_name,
        logger=logger,
        max_workers=max_workers,
    )

    llm_selected_df = ranked_df[ranked_df["LLM_SCORE"] >= 80].copy()
    if llm_selected_df.empty:
        logger.warning("No rows scored >= 80 — keeping top-ranked rows.")
        llm_selected_df = ranked_df.head(
            min(20, len(ranked_df))
        ).copy()

    logger.info(
        "Rows after LLM ranking: %s", f"{len(llm_selected_df):,}"
    )
    return llm_selected_df, ranked_df, disease_variants
