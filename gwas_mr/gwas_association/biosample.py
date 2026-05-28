"""Biosample / tissue normalisation, matching, and LLM-based filtering."""

import json
import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Tuple

import pandas as pd
from openai import OpenAI

from .config import MODEL_NAME, STEP6_PROGRESS_EVERY, TEXTUAL_BIOSAMPLE_HINTS
from .helpers import normalize_text
from .llm_utils import llm_text


def keyword_biosample_match(text: str, biosample: str) -> bool:
    """Check whether *text* mentions *biosample* or a known synonym.

    Args:
        text: Free-text field to search (e.g. study description).
        biosample: Canonical biosample term.

    Returns:
        *True* if a direct or synonym match is found.
    """
    text_n = normalize_text(text)
    bio_n = normalize_text(biosample)

    if bio_n in text_n:
        return True

    for canonical, hints in TEXTUAL_BIOSAMPLE_HINTS.items():
        if bio_n == canonical or bio_n in hints:
            for hint in hints:
                if hint in text_n:
                    return True
    return False


def normalize_biosample_input(
    biosample: str,
    client: OpenAI,
    model_name: str = MODEL_NAME,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Map a user-provided biosample term to a canonical label.

    First checks the static
    :data:`~gwas_association.config.TEXTUAL_BIOSAMPLE_HINTS` table.  Falls back
    to an LLM prompt for unrecognised terms.

    Args:
        biosample: Raw user biosample / tissue string.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        logger: Optional logger instance.

    Returns:
        Canonical biosample label (lowercase).
    """
    biosample_n = normalize_text(biosample)
    if not biosample_n:
        return biosample

    for canonical, hints in TEXTUAL_BIOSAMPLE_HINTS.items():
        all_terms = [canonical] + hints
        if biosample_n in [normalize_text(term) for term in all_terms]:
            return canonical

    prompt = f"""
You are a biomedical tissue and biosample normalizer.

Task:
Map the user's biosample or tissue term to the closest canonical label.

Allowed canonical labels:
- pbmc
- whole blood
- blood
- biopsy
- tissue
- skin
- kidney
- lung
- serum
- plasma

If none fit well, return the original term normalized to lowercase.

Return ONLY valid JSON:
{{
  "canonical_biosample": "..."
}}

User biosample: {biosample}
""".strip()

    try:
        raw = llm_text(
            prompt, client=client, model_name=model_name,
            max_output_tokens=80,
        )
        parsed = json.loads(raw)
        canonical = str(parsed.get("canonical_biosample", "")).strip()
        return canonical if canonical else biosample_n
    except Exception:
        if logger:
            logger.warning(
                "Biosample normalisation via LLM failed; using raw term."
            )
        return biosample_n


def simple_tissue_filter(df: pd.DataFrame, tissue: str) -> pd.DataFrame:
    """Keep rows where *tissue* appears in the study or sample-size columns.

    Args:
        df: Pre-filtered GWAS DataFrame.
        tissue: Normalised tissue term.

    Returns:
        Filtered copy.
    """
    tissue_n = normalize_text(tissue)
    if not tissue_n:
        return df.copy()

    mask = (
        df["INITIAL SAMPLE SIZE"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(re.escape(tissue_n), regex=True)
        | df["STUDY"]
        .fillna("")
        .astype(str)
        .str.lower()
        .str.contains(re.escape(tissue_n), regex=True)
    )

    return df[mask].copy()


def enforce_disease_tissue_consistency(
    df: pd.DataFrame,
    disease: str,
    tissue: str,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Reject obvious disease-tissue mismatches via hard-coded rules.

    For example, a *cutaneous* disease requires *skin* tissue.

    Args:
        df: Pre-filtered GWAS DataFrame.
        disease: Corrected disease name.
        tissue: Normalised tissue term.
        logger: Optional logger instance.

    Returns:
        *df* unchanged, or an empty DataFrame if a mismatch is detected.
    """
    disease_n = normalize_text(disease)
    tissue_n = normalize_text(tissue)

    rules = [
        ("cutaneous", "skin"),
        ("skin", "skin"),
        ("lung", "lung"),
        ("pulmonary", "lung"),
        ("oral", "mouth"),
        ("head", "mouth"),
        ("head and neck", "mouth"),
        ("esophageal", "esophagus"),
    ]

    for keyword, required_tissue in rules:
        if keyword in disease_n and required_tissue not in tissue_n:
            if logger:
                logger.warning(
                    "Mismatch: disease implies '%s' but user selected '%s'",
                    required_tissue, tissue,
                )
            return df.iloc[0:0].copy()

    return df


def _check_tissue_relevance(
    pos: int,
    row: pd.Series,
    disease: str,
    tissue: str,
    client: OpenAI,
    model_name: str,
) -> Tuple[int, Optional[pd.Series]]:
    """Check one row for tissue relevance (runs in a worker thread)."""
    text = (
        str(row.get("STUDY", ""))
        + " "
        + str(row.get("INITIAL SAMPLE SIZE", ""))
    )

    prompt = f"""
You are a biomedical expert.

Disease: {disease}
Tissue: {tissue}

Study text:
{text}

Important:
- If disease subtype implies a different tissue (for example cutaneous = skin, lung/pulmonary = lung), reject it.
- Only accept if BOTH disease subtype AND tissue match.

Question:
Is this study relevant to BOTH disease AND tissue?

Answer ONLY JSON:
{{"relevant": true/false}}
""".strip()

    try:
        raw = llm_text(
            prompt, client=client, model_name=model_name,
            max_output_tokens=50,
        )
        result = json.loads(raw)
        if result.get("relevant") is True:
            return pos, row.copy()
        return pos, None
    except Exception:
        return pos, None


def llm_filter_tissue(
    df: pd.DataFrame,
    disease: str,
    tissue: str,
    client: OpenAI,
    model_name: str = MODEL_NAME,
    logger: Optional[logging.Logger] = None,
    progress_every: int = STEP6_PROGRESS_EVERY,
    max_workers: int = 10,
) -> pd.DataFrame:
    """Use an LLM to keep only rows relevant to both *disease* and *tissue*.

    Rows are checked **concurrently** using a thread pool (controlled by
    *max_workers*) for significantly faster execution.

    Args:
        df: Pre-filtered GWAS DataFrame.
        disease: Corrected disease name.
        tissue: Normalised tissue term.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        logger: Optional logger instance.
        progress_every: Log progress every *N* completed rows.
        max_workers: Number of concurrent LLM threads.

    Returns:
        DataFrame containing only LLM-approved rows.
    """
    total_rows = len(df)

    if logger:
        logger.info(
            "Checking tissue relevance for %d GWAS rows "
            "(%d concurrent workers)...",
            total_rows, max_workers,
        )

    work_items = [
        (pos, row)
        for pos, (_, row) in enumerate(df.iterrows(), start=1)
    ]

    selected_rows = []
    completed = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _check_tissue_relevance,
                pos, row, disease, tissue, client, model_name,
            ): pos
            for pos, row in work_items
        }

        for future in as_completed(futures):
            pos, kept_row = future.result()
            completed += 1

            if kept_row is not None:
                selected_rows.append(kept_row)

            if logger and (
                completed == 1
                or completed % progress_every == 0
                or completed == total_rows
            ):
                logger.info(
                    "  Tissue filter: %d/%d checked, %d kept",
                    completed, total_rows, len(selected_rows),
                )

    if not selected_rows:
        return pd.DataFrame(columns=df.columns)

    return pd.DataFrame(selected_rows)
