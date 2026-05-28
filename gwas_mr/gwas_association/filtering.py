"""GWAS disease / trait filtering utilities.

Provides strict keyword filters, similarity-based preview / ranking, and
RAG-index-accelerated partial matching.
"""

import difflib
import re
from typing import Dict, List, Optional

import pandas as pd

from .helpers import get_trait_match_column, normalize_text, tokenize_for_index


def get_all_disease_values(gwas_df: pd.DataFrame) -> List[str]:
    """Return every non-null trait string from the GWAS trait column.

    Args:
        gwas_df: Full GWAS association DataFrame.

    Returns:
        List of raw trait strings (stripped, not normalised).
    """
    trait_col = get_trait_match_column(gwas_df)
    return (
        gwas_df[trait_col]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )


def preview_closest_traits(
    gwas_df: pd.DataFrame,
    disease_input: str,
    rag_index: Optional[Dict[str, object]] = None,
    top_n: int = 15,
) -> pd.DataFrame:
    """Score and rank GWAS traits by similarity to *disease_input*.

    Uses a combination of sequence similarity, token overlap, and substring
    boosting.  An optional RAG index narrows the candidate set.

    Args:
        gwas_df: Full GWAS association DataFrame.
        disease_input: User-provided disease name.
        rag_index: Pre-built RAG index (speeds up candidate selection).
        top_n: Maximum number of preview rows to return.

    Returns:
        DataFrame with columns ``TRAIT_COLUMN_USED``, ``NORMALIZED_TRAIT``,
        ``DISEASE_TRAIT``, ``SIMILARITY_SCORE``, and ``ROW_COUNT``, sorted
        by descending score.
    """
    disease_n = normalize_text(disease_input)
    trait_col = get_trait_match_column(gwas_df)

    if rag_index:
        cached_trait_col = rag_index.get("trait_column")
        if cached_trait_col and cached_trait_col in gwas_df.columns:
            trait_col = cached_trait_col

    trait_series = gwas_df[trait_col].dropna().astype(str).str.strip()
    normalized_to_originals: Dict[str, List[str]] = {}
    normalized_counts: Dict[str, int] = {}

    for trait in trait_series.tolist():
        trait_n = normalize_text(trait)
        if not trait_n:
            continue
        normalized_to_originals.setdefault(trait_n, []).append(trait)
        normalized_counts[trait_n] = normalized_counts.get(trait_n, 0) + 1

    if rag_index and rag_index.get("unique_traits"):
        candidate_traits = rag_index.get("unique_traits", [])
    else:
        candidate_traits = list(normalized_to_originals.keys())

    scored_rows = []
    for trait in candidate_traits:
        trait_n = normalize_text(trait)
        similarity = difflib.SequenceMatcher(None, disease_n, trait_n).ratio()
        token_overlap = len(
            set(tokenize_for_index(disease_n))
            & set(tokenize_for_index(trait_n))
        )
        boost = 0.0

        if disease_n and disease_n in trait_n:
            boost += 0.35
        if token_overlap:
            boost += min(0.25, token_overlap * 0.08)

        final_score = min(1.0, similarity + boost)
        original_labels = normalized_to_originals.get(
            trait_n, [str(trait).strip()]
        )
        label_counts = pd.Series(original_labels).value_counts()
        trait_display = str(label_counts.index[0]).strip()
        row_count = int(normalized_counts.get(trait_n, 0))

        scored_rows.append({
            "TRAIT_COLUMN_USED": trait_col,
            "NORMALIZED_TRAIT": trait_n,
            "DISEASE_TRAIT": trait_display,
            "SIMILARITY_SCORE": round(final_score, 4),
            "ROW_COUNT": row_count,
        })

    preview_df = pd.DataFrame(scored_rows)
    if preview_df.empty:
        return preview_df

    preview_df = preview_df.sort_values(
        by=["SIMILARITY_SCORE", "ROW_COUNT", "DISEASE_TRAIT"],
        ascending=[False, False, True],
    ).head(top_n).reset_index(drop=True)

    return preview_df


def choose_best_preview_traits(
    preview_df: pd.DataFrame,
    min_score: float = 0.80,
    max_traits: int = 8,
) -> List[str]:
    """Select the best normalised traits from a preview DataFrame.

    Args:
        preview_df: Output of :func:`preview_closest_traits`.
        min_score: Minimum similarity score threshold.
        max_traits: Maximum traits to return.

    Returns:
        List of normalised trait strings.
    """
    if preview_df.empty:
        return []

    chosen = preview_df[
        preview_df["SIMILARITY_SCORE"] >= min_score
    ].head(max_traits)
    if chosen.empty:
        chosen = preview_df.head(min(3, len(preview_df)))

    return chosen["NORMALIZED_TRAIT"].astype(str).tolist()


def filter_rows_by_preview_traits(
    gwas_df: pd.DataFrame,
    normalized_traits: List[str],
    trait_col: str,
) -> pd.DataFrame:
    """Keep only rows whose trait column matches one of *normalized_traits*.

    Args:
        gwas_df: GWAS DataFrame to filter.
        normalized_traits: Allowed normalised trait values.
        trait_col: Column name to match against.

    Returns:
        Filtered copy of *gwas_df*.
    """
    if not normalized_traits:
        return gwas_df.copy()

    allowed = {normalize_text(x) for x in normalized_traits if str(x).strip()}
    mask = gwas_df[trait_col].fillna("").astype(str).apply(
        lambda x: normalize_text(x) in allowed
    )
    return gwas_df[mask].copy()


def strict_disease_filter(
    gwas_df: pd.DataFrame,
    disease: str,
    trait_col: str,
) -> pd.DataFrame:
    """Filter GWAS rows by strict keyword match, excluding mixed studies.

    Rows whose trait column contains ``mtag``, ``combined``, or the
    conjunctions ``or`` / ``and`` (as whole words) are excluded to avoid
    multi-phenotype studies.

    Args:
        gwas_df: Full GWAS DataFrame.
        disease: Disease search term.
        trait_col: Column name to match against.

    Returns:
        Filtered copy of *gwas_df*.
    """
    disease_n = normalize_text(disease)

    if "aml" in disease_n or "acute myeloid leukemia" in disease_n:
        keywords = [
            "acute myeloid leukemia",
            "aml",
            "core binding factor acute myeloid leukemia",
            "acute myelogenous leukemia",
        ]
    else:
        keywords = [disease_n]

    mask = (
        gwas_df[trait_col]
        .fillna("")
        .astype(str)
        .str.lower()
        .apply(lambda x: any(k in x for k in keywords))
    )

    df = gwas_df[mask].copy()
    df = df[
        ~df[trait_col]
        .fillna("")
        .astype(str)
        .str.contains(
            r"mtag|combined|\bor\b|\band\b",
            case=False,
            regex=True,
            na=False,
        )
    ]

    return df


def simple_trait_filter(
    gwas_df: pd.DataFrame,
    user_query: str,
    trait_col: str,
) -> pd.DataFrame:
    """Filter GWAS rows by simple substring match on the trait column.

    Args:
        gwas_df: Full GWAS DataFrame.
        user_query: Disease term to search for.
        trait_col: Column to search in.

    Returns:
        Matching rows (copy).
    """
    query = normalize_text(user_query)
    if not query:
        return gwas_df.iloc[0:0].copy()

    return gwas_df[
        gwas_df[trait_col]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(re.escape(query), regex=True)
    ].copy()


def exact_and_partial_match(
    gwas_df: pd.DataFrame,
    disease_terms: List[str],
    rag_index: Optional[Dict[str, object]] = None,
) -> pd.DataFrame:
    """Combine exact and partial substring matches across disease terms.

    When a RAG index is available, candidate traits are resolved via token
    lookup first for performance.

    Args:
        gwas_df: Full GWAS DataFrame.
        disease_terms: One or more disease terms to match.
        rag_index: Pre-built RAG index.

    Returns:
        Concatenated, de-duplicated DataFrame with an additional
        ``MATCHED_DISEASE_TERM`` column.
    """
    disease_terms_n = [
        normalize_text(x) for x in disease_terms if str(x).strip()
    ]
    out_frames = []

    for term in disease_terms_n:
        matched_indices: set = set()

        if rag_index:
            token_to_traits = rag_index.get("token_to_traits", {})
            trait_to_rows = rag_index.get("trait_to_rows", {})

            term_tokens = tokenize_for_index(term)
            candidate_traits: set = set()

            for token in term_tokens:
                candidate_traits.update(token_to_traits.get(token, []))

            if not candidate_traits:
                candidate_traits = set(rag_index.get("unique_traits", []))

            for trait in candidate_traits:
                if term in trait:
                    matched_indices.update(trait_to_rows.get(trait, []))

        if matched_indices:
            temp = gwas_df.iloc[sorted(matched_indices)].copy()
        else:
            mask = (
                gwas_df["DISEASE/TRAIT"]
                .fillna("")
                .astype(str)
                .str.lower()
                .str.contains(re.escape(term), regex=True)
            )
            temp = gwas_df[mask].copy()

        if not temp.empty:
            temp["MATCHED_DISEASE_TERM"] = term
            out_frames.append(temp)

    if not out_frames:
        return pd.DataFrame(
            columns=list(gwas_df.columns) + ["MATCHED_DISEASE_TERM"]
        )

    merged = pd.concat(out_frames, ignore_index=True).drop_duplicates()
    return merged


def fallback_rank_candidates(
    gwas_df: pd.DataFrame,
    disease_input: str,
    top_n: int = 200,
) -> pd.DataFrame:
    """Last-resort fuzzy ranking when strict / partial filters return nothing.

    Args:
        gwas_df: Full GWAS DataFrame.
        disease_input: User-provided disease name.
        top_n: Maximum rows to return.

    Returns:
        Up to *top_n* rows from the best-matching traits.
    """
    trait_col = get_trait_match_column(gwas_df)
    unique_traits = gwas_df[trait_col].dropna().astype(str).unique().tolist()
    scored = []

    for trait in unique_traits:
        score = difflib.SequenceMatcher(
            None, normalize_text(disease_input), normalize_text(trait)
        ).ratio()
        scored.append((trait, score))

    scored = sorted(scored, key=lambda x: x[1], reverse=True)[:25]
    top_traits = {trait for trait, score in scored if score >= 0.35}

    temp = gwas_df[gwas_df[trait_col].astype(str).isin(top_traits)].copy()
    return temp.head(top_n)
