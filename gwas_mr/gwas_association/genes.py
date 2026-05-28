"""Gene column detection, token parsing, and GWAS gene-overlap mapping."""

import re
from typing import List, Optional, Set

import pandas as pd

from .config import GENE_COLUMN_CANDIDATES
from .helpers import normalize_text


def find_gene_column(
    df: pd.DataFrame,
    gene_column_candidates: Optional[List[str]] = None,
) -> str:
    """Identify the gene-symbol column in *df*.

    Tries an exact normalised match against the candidate list first, then
    falls back to heuristic pattern matching (column name contains ``gene``
    and ``symbol``).

    Args:
        df: Input transcriptomics DataFrame.
        gene_column_candidates: Column names to try.  Falls back to
            :data:`~gwas_association.config.GENE_COLUMN_CANDIDATES`.

    Returns:
        The matched column name as it appears in *df*.

    Raises:
        ValueError: If no gene column can be identified.
    """
    candidates = gene_column_candidates or GENE_COLUMN_CANDIDATES
    normalized = {
        col: normalize_text(col).replace("_", " ") for col in df.columns
    }

    candidate_norms = [normalize_text(c).replace("_", " ") for c in candidates]
    for col, norm_col in normalized.items():
        if norm_col in candidate_norms:
            return col

    for col, norm_col in normalized.items():
        if "gene" in norm_col and ("symbol" in norm_col or norm_col == "gene"):
            return col

    raise ValueError(
        f"No gene column found.  Columns available: {list(df.columns)}"
    )


def parse_gene_tokens(raw_gene_value: str) -> List[str]:
    """Split a raw gene-field value into individual uppercase gene symbols.

    Delimiters recognised: ``,``, ``;``, ``/``, ``|``, newline.
    Placeholder values (*NR*, *N/A*, *NA*, ``-``, ``--``) are dropped.

    Args:
        raw_gene_value: Raw cell value (may be NaN).

    Returns:
        De-duplicated list of uppercase gene tokens.
    """
    if pd.isna(raw_gene_value):
        return []

    text = str(raw_gene_value).strip()
    if not text:
        return []

    parts = re.split(r"[,;/\|\n]+", text)
    cleaned: List[str] = []

    for p in parts:
        token = p.strip().upper()
        token = re.sub(r"\s+", " ", token)
        if token and token not in {"NR", "N/A", "NA", "-", "--"}:
            cleaned.append(token)

    return list(dict.fromkeys(cleaned))


def build_gene_set_from_input(df: pd.DataFrame, gene_col: str) -> Set[str]:
    """Collect all unique gene symbols from *gene_col* in *df*.

    Args:
        df: Input transcriptomics DataFrame.
        gene_col: Name of the column containing gene symbols.

    Returns:
        Set of uppercase gene symbols.
    """
    genes: Set[str] = set()
    for value in df[gene_col].dropna():
        for g in parse_gene_tokens(value):
            genes.add(g)
    return genes


def extract_all_possible_genes(row: pd.Series) -> List[str]:
    """Extract gene symbols from every gene-related GWAS field in *row*.

    Fields inspected: ``REPORTED GENE(S)``, ``MAPPED_GENE``,
    ``UPSTREAM_GENE_ID``, ``DOWNSTREAM_GENE_ID``, ``SNP_GENE_IDS``.

    Args:
        row: A single GWAS association row.

    Returns:
        De-duplicated list of uppercase gene symbols.
    """
    gene_fields = [
        "REPORTED GENE(S)",
        "MAPPED_GENE",
        "UPSTREAM_GENE_ID",
        "DOWNSTREAM_GENE_ID",
        "SNP_GENE_IDS",
    ]

    all_genes: List[str] = []
    for field in gene_fields:
        value = row.get(field, "")
        all_genes.extend(parse_gene_tokens(value))

    return list(dict.fromkeys(all_genes))


def gene_overlap(
    reported_genes: List[str],
    input_genes: Set[str],
) -> List[str]:
    """Return genes present in both *reported_genes* and *input_genes*.

    Comparison is case-insensitive; original case from *reported_genes* is
    preserved in the output.

    Args:
        reported_genes: Genes extracted from a GWAS row.
        input_genes: Set of patient / input gene symbols.

    Returns:
        De-duplicated list of overlapping genes.
    """
    input_genes_upper = {g.upper() for g in input_genes}
    matches = [g for g in reported_genes if g.upper() in input_genes_upper]
    return list(dict.fromkeys(matches))


def map_input_genes_to_gwas(
    filtered_df: pd.DataFrame,
    input_genes: Set[str],
) -> pd.DataFrame:
    """Map patient genes against GWAS rows, keeping only overlapping rows.

    Two columns are appended to each matched row:

    * ``INPUT_GENE_MATCHES`` -- semicolon-separated matched gene symbols.
    * ``NUM_INPUT_GENE_MATCHES`` -- count of matched genes.

    Args:
        filtered_df: Pre-filtered GWAS DataFrame.
        input_genes: Set of patient / input gene symbols.

    Returns:
        DataFrame containing only rows with at least one gene overlap.
    """
    mapped_rows = []

    for _, row in filtered_df.iterrows():
        row_genes = extract_all_possible_genes(row)
        overlap = gene_overlap(row_genes, input_genes)

        if overlap:
            row_copy = row.copy()
            row_copy["INPUT_GENE_MATCHES"] = "; ".join(overlap)
            row_copy["NUM_INPUT_GENE_MATCHES"] = len(overlap)
            mapped_rows.append(row_copy)

    if not mapped_rows:
        return pd.DataFrame(
            columns=list(filtered_df.columns)
            + ["INPUT_GENE_MATCHES", "NUM_INPUT_GENE_MATCHES"]
        )

    return pd.DataFrame(mapped_rows)


def expand_mapped_rows_per_gene(mapped_df: pd.DataFrame) -> pd.DataFrame:
    """Explode *mapped_df* so that each matched gene gets its own row.

    A new column ``MAPPED_INPUT_GENE`` holds the individual gene symbol.

    Args:
        mapped_df: Output of :func:`map_input_genes_to_gwas`.

    Returns:
        Expanded DataFrame with one row per gene match.
    """
    expanded = []

    for _, row in mapped_df.iterrows():
        matches = str(row.get("INPUT_GENE_MATCHES", "")).split(";")
        for gene in matches:
            gene = gene.strip()
            if gene:
                row_copy = row.copy()
                row_copy["MAPPED_INPUT_GENE"] = gene
                expanded.append(row_copy)

    if not expanded:
        return pd.DataFrame(
            columns=list(mapped_df.columns) + ["MAPPED_INPUT_GENE"]
        )

    return pd.DataFrame(expanded)
