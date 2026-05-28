"""Output file writers, summary generators, and matplotlib-based plots."""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd

from .genes import extract_all_possible_genes


def save_summary_json(
    output_path: Path,
    user_disease: str,
    corrected_disease: str,
    biosample: str,
    disease_variants: List[str],
    index_rows_total: int,
    disease_candidate_rows: int,
    llm_selected_rows: int,
    mapped_gene_rows: int,
    unique_mapped_genes: int,
    input_file: str,
    gene_column: str,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Write a machine-readable JSON run summary.

    Args:
        output_path: Destination file path.
        user_disease: Original user input.
        corrected_disease: Spelling-corrected disease name.
        biosample: Normalised biosample term.
        disease_variants: Variant terms used during filtering.
        index_rows_total: Total rows in the GWAS index.
        disease_candidate_rows: Rows after disease filtering.
        llm_selected_rows: Rows after LLM selection.
        mapped_gene_rows: Rows with gene overlaps.
        unique_mapped_genes: Count of unique matched genes.
        input_file: Name of the input transcriptomics file.
        gene_column: Auto-detected gene column name.
        logger: Optional logger instance.
    """
    summary = {
        "user_disease": user_disease,
        "corrected_disease": corrected_disease,
        "biosample": biosample,
        "disease_variants_used": disease_variants,
        "index_rows_total": index_rows_total,
        "disease_candidate_rows": disease_candidate_rows,
        "llm_selected_rows": llm_selected_rows,
        "mapped_gene_rows": mapped_gene_rows,
        "unique_mapped_genes": unique_mapped_genes,
        "input_file": input_file,
        "gene_column_detected": gene_column,
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    if logger:
        logger.info("Saved summary JSON: %s", output_path.name)


def save_summary_txt(
    output_path: Path,
    user_disease: str,
    corrected_disease: str,
    biosample: str,
    disease_candidate_rows: int,
    llm_selected_rows: int,
    mapped_gene_rows: int,
    unique_mapped_genes: int,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Write a human-readable plain-text run summary.

    Args:
        output_path: Destination file path.
        user_disease: Original user input.
        corrected_disease: Spelling-corrected disease name.
        biosample: Normalised biosample term.
        disease_candidate_rows: Rows after disease filtering.
        llm_selected_rows: Rows after LLM selection.
        mapped_gene_rows: Rows with gene overlaps.
        unique_mapped_genes: Count of unique matched genes.
        logger: Optional logger instance.
    """
    lines = [
        f"User disease: {user_disease}",
        f"Corrected disease: {corrected_disease}",
        f"Tissue/Biosample: {biosample}",
        f"Total GWAS candidate rows: {disease_candidate_rows}",
        f"Selected GWAS rows after filtering: {llm_selected_rows}",
        f"Mapped gene rows: {mapped_gene_rows}",
        f"Unique mapped genes: {unique_mapped_genes}",
    ]

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    if logger:
        logger.info("Saved summary TXT: %s", output_path.name)


def save_llm_reasoning(
    df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Write the top-50 LLM-scored rows to a CSV for auditability.

    Args:
        df: DataFrame expected to contain ``LLM_SCORE`` and ``LLM_REASON``.
        output_path: Destination CSV path.
        logger: Optional logger instance.
    """
    cols = ["DISEASE/TRAIT", "STUDY", "LLM_SCORE", "LLM_REASON"]
    available_cols = [c for c in cols if c in df.columns]

    if not available_cols:
        pd.DataFrame().to_csv(output_path, index=False)
        return

    df[available_cols].head(50).to_csv(output_path, index=False)

    if logger:
        logger.info("Saved LLM reasoning: %s", output_path.name)


def save_run_overview(
    output_path: Path,
    user_disease: str,
    corrected_disease: str,
    biosample: str,
    input_file: str,
    main_result_file: Path,
    reasoning_file: Path,
    transcriptomics_file: Path,
    expanded_file: Path,
    candidate_file: Path,
    summary_json_file: Path,
    summary_txt_file: Path,
    gene_plot_file: Path,
    trait_plot_file: Path,
    population_plot_file: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Write a run-overview text file listing all output artefacts.

    Args:
        output_path: Destination file path.
        user_disease: Original user input.
        corrected_disease: Spelling-corrected disease name.
        biosample: Normalised biosample term.
        input_file: Name of the input transcriptomics file.
        main_result_file: Primary GWAS output CSV.
        reasoning_file: LLM reasoning CSV.
        transcriptomics_file: Gene-overlap CSV.
        expanded_file: Per-gene expanded CSV.
        candidate_file: Pre-ranking candidate CSV.
        summary_json_file: JSON summary.
        summary_txt_file: Text summary.
        gene_plot_file: Gene bar-chart PNG.
        trait_plot_file: Trait bar-chart PNG.
        population_plot_file: Population bar-chart PNG.
        logger: Optional logger instance.
    """
    lines = [
        "RUN OVERVIEW",
        "",
        f"Input disease: {user_disease}",
        f"Corrected disease used: {corrected_disease}",
        f"Biosample/Tissue: {biosample}",
        f"Input transcriptomics file used: {input_file}",
        "",
        "MAIN FILE TO OPEN FIRST:",
        f"{main_result_file.name}",
        "This is the final GWAS result chosen by the pipeline.",
        "",
        "OTHER IMPORTANT FILES:",
        f"{reasoning_file.name} -> why rows were ranked and selected",
        f"{transcriptomics_file.name} -> final GWAS rows that overlap "
        "with your input genes",
        f"{expanded_file.name} -> one row per matched gene",
        "",
        "SUPPORTING FILES:",
        f"{candidate_file.name} -> disease candidates before LLM ranking",
        f"{summary_json_file.name} -> machine-readable summary",
        f"{summary_txt_file.name} -> human-readable summary",
        f"{gene_plot_file.name} -> top mapped genes plot",
        f"{trait_plot_file.name} -> selected traits plot",
        f"{population_plot_file.name} -> population context plot",
    ]

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    if logger:
        logger.info("Saved run overview: %s", output_path.name)


def make_knowledge_graph(
    disease_name: str,
    selected_df: pd.DataFrame,
    mapped_expanded_df: pd.DataFrame,
    output_path: Path,
    max_genes: int = 30,
    max_snps: int = 60,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Render a Disease -> Gene -> SNP knowledge graph as an interactive HTML.

    The output is a standalone HTML file using vis.js (via pyvis) that
    supports zoom, pan, and drag so overlapping nodes can be separated
    by the viewer.

    Args:
        disease_name: Central disease node label.
        selected_df: GWAS rows used for SNP extraction.
        mapped_expanded_df: Per-gene expanded DataFrame for gene counts.
        output_path: Destination file path.  If it ends with ``.png``,
            the extension is replaced with ``.html`` automatically.
        max_genes: Maximum gene nodes to display.
        max_snps: Maximum SNP nodes to display.
        logger: Optional logger instance.
    """
    from pyvis.network import Network

    if selected_df.empty:
        if logger:
            logger.info("Skipping knowledge graph — no selected rows.")
        return

    output_path = Path(output_path)
    if output_path.suffix.lower() == ".png":
        output_path = output_path.with_suffix(".html")

    gene_counts: Dict[str, int] = {}
    if (
        not mapped_expanded_df.empty
        and "MAPPED_INPUT_GENE" in mapped_expanded_df.columns
    ):
        gene_counts = (
            mapped_expanded_df["MAPPED_INPUT_GENE"]
            .value_counts()
            .head(max_genes)
            .to_dict()
        )

    genes = list(gene_counts.keys())
    snp_links: List[tuple] = []
    used_snps: set = set()

    for _, row in selected_df.iterrows():
        snp = str(row.get("SNPS", "")).strip()
        if not snp or snp.lower() == "nan" or snp in used_snps:
            continue

        linked_gene = None
        row_genes = extract_all_possible_genes(row)
        for gene in row_genes:
            if gene in genes:
                linked_gene = gene
                break

        snp_links.append((snp, linked_gene))
        used_snps.add(snp)
        if len(snp_links) >= max_snps:
            break

    net = Network(
        height="900px",
        width="100%",
        bgcolor="#ffffff",
        font_color="#333333",
        directed=False,
        cdn_resources="remote",
    )

    net.barnes_hut(
        gravity=-3000,
        central_gravity=0.3,
        spring_length=150,
        spring_strength=0.01,
        damping=0.09,
    )

    net.add_node(
        disease_name,
        label=disease_name,
        color="#d73027",
        size=45,
        font={"size": 16, "color": "white"},
        shape="dot",
        title=f"Disease: {disease_name}",
    )

    for gene, count in gene_counts.items():
        size = 20 + min(count * 3, 25)
        net.add_node(
            gene,
            label=gene,
            color="#74add1",
            size=size,
            font={"size": 13, "color": "#000000"},
            shape="dot",
            title=f"Gene: {gene}\nGWAS associations: {count}",
        )
        net.add_edge(
            disease_name, gene,
            color="#4575b4",
            width=2,
        )

    for snp, linked_gene in snp_links:
        net.add_node(
            snp,
            label=snp,
            color="#fee090",
            size=12,
            font={"size": 10, "color": "#333333"},
            shape="dot",
            borderWidth=1,
            title=f"SNP: {snp}",
        )
        if linked_gene and linked_gene in gene_counts:
            net.add_edge(
                linked_gene, snp,
                color="#fdae61",
                width=1,
            )
        else:
            net.add_edge(
                disease_name, snp,
                color="#fdae61",
                width=1,
                dashes=True,
            )

    net.set_options("""
    {
      "interaction": {
        "hover": true,
        "navigationButtons": true,
        "keyboard": true,
        "zoomView": true,
        "dragNodes": true
      },
      "physics": {
        "enabled": true,
        "barnesHut": {
          "gravitationalConstant": -3000,
          "centralGravity": 0.3,
          "springLength": 150,
          "springConstant": 0.01,
          "damping": 0.09,
          "avoidOverlap": 0.5
        },
        "stabilization": {
          "enabled": true,
          "iterations": 200
        }
      }
    }
    """)

    net.save_graph(str(output_path))

    import re as _re
    import shutil as _shutil

    html_text = output_path.read_text(encoding="utf-8")
    html_text = _re.sub(
        r'<script\s+src="lib/bindings/utils\.js"\s*></script>\s*', "", html_text
    )
    html_text = _re.sub(
        r"<center>\s*<h1>.*?</h1>\s*</center>\s*", "", html_text
    )
    title_html = (
        '<center><h1 style="margin:20px 0">Disease-Gene-SNP Knowledge Graph: '
        f'{disease_name}</h1></center>\n'
    )
    html_text = html_text.replace("<body>", f"<body>\n{title_html}", 1)
    output_path.write_text(html_text, encoding="utf-8")

    lib_dir = output_path.parent / "lib"
    if lib_dir.is_dir():
        _shutil.rmtree(lib_dir, ignore_errors=True)

    if logger:
        logger.info("Saved interactive knowledge graph: %s", output_path.name)


def make_gene_plot(
    mapped_expanded_df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Bar chart of the top-20 most frequently matched GWAS genes.

    Args:
        mapped_expanded_df: Per-gene expanded DataFrame.
        output_path: Destination PNG path.
        logger: Optional logger instance.
    """
    if (
        mapped_expanded_df.empty
        or "MAPPED_INPUT_GENE" not in mapped_expanded_df.columns
    ):
        return

    top_counts = (
        mapped_expanded_df["MAPPED_INPUT_GENE"].value_counts().head(20)
    )
    if top_counts.empty:
        return

    plt.figure(figsize=(12, 6))
    top_counts.plot(kind="bar")
    plt.title("Top Mapped GWAS Genes")
    plt.xlabel("Gene")
    plt.ylabel("Number of Matching GWAS Rows")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    if logger:
        logger.info("Saved gene plot: %s", output_path.name)


def make_trait_plot(
    filtered_df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Bar chart of the top-15 disease / trait categories in the filtered set.

    Args:
        filtered_df: GWAS DataFrame after filtering.
        output_path: Destination PNG path.
        logger: Optional logger instance.
    """
    if filtered_df.empty:
        return

    top_traits = (
        filtered_df["DISEASE/TRAIT"].astype(str).value_counts().head(15)
    )
    if top_traits.empty:
        return

    plt.figure(figsize=(12, 6))
    top_traits.plot(kind="bar")
    plt.title("Top Selected GWAS Disease/Trait Categories")
    plt.xlabel("Disease/Trait")
    plt.ylabel("Number of Rows")
    plt.xticks(rotation=60, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    if logger:
        logger.info("Saved trait plot: %s", output_path.name)


def make_population_plot(
    filtered_df: pd.DataFrame,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Bar chart of GWAS rows grouped by population ancestry keywords.

    Args:
        filtered_df: GWAS DataFrame after filtering.
        output_path: Destination PNG path.
        logger: Optional logger instance.
    """
    if (
        filtered_df.empty
        or "INITIAL SAMPLE SIZE" not in filtered_df.columns
    ):
        return

    ancestry_keywords = {
        "European": 0,
        "Asian": 0,
        "African": 0,
        "Hispanic": 0,
        "Mixed/Other": 0,
    }

    for val in filtered_df["INITIAL SAMPLE SIZE"].fillna("").astype(str):
        val_low = val.lower()
        if "european" in val_low:
            ancestry_keywords["European"] += 1
        elif "asian" in val_low:
            ancestry_keywords["Asian"] += 1
        elif "african" in val_low:
            ancestry_keywords["African"] += 1
        elif "hispanic" in val_low or "latin" in val_low:
            ancestry_keywords["Hispanic"] += 1
        else:
            ancestry_keywords["Mixed/Other"] += 1

    series = pd.Series(ancestry_keywords)
    if series.sum() == 0:
        return

    plt.figure(figsize=(8, 5))
    series.plot(kind="bar")
    plt.title("GWAS Association by Population Context")
    plt.xlabel("Population Group")
    plt.ylabel("Row Count")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    if logger:
        logger.info("Saved population plot: %s", output_path.name)
