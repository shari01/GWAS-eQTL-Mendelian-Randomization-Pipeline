"""LLM-powered helpers: disease correction, variant generation, row ranking."""

import difflib
import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple

import pandas as pd
from openai import OpenAI

from .config import MODEL_NAME, STEP6_PROGRESS_EVERY
from .helpers import normalize_text


def llm_text(
    prompt: str,
    client: OpenAI,
    model_name: str = MODEL_NAME,
    max_output_tokens: int = 500,
) -> str:
    """Send a single prompt to the OpenAI API and return the trimmed text.

    Args:
        prompt: The prompt string.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        max_output_tokens: Response length cap.

    Returns:
        Stripped response text.
    """
    response = client.responses.create(
        model=model_name,
        input=prompt,
        max_output_tokens=max_output_tokens,
    )
    return response.output_text.strip()


def correct_disease_name(
    user_disease: str,
    all_disease_values: List[str],
    client: OpenAI,
    model_name: str = MODEL_NAME,
    logger: Optional[logging.Logger] = None,
) -> str:
    """Correct potential spelling errors in a user-provided disease name.

    First attempts a fast fuzzy match against known trait values; falls back
    to an LLM prompt when no close match is found.

    Args:
        user_disease: Raw disease name from the user.
        all_disease_values: Known trait strings from the GWAS index.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        logger: Optional logger instance.

    Returns:
        Corrected disease name (unchanged if correction fails).
    """
    user_disease_n = normalize_text(user_disease)
    disease_lookup = list(dict.fromkeys(
        [str(x) for x in all_disease_values if str(x).strip()]
    ))

    close_matches = difflib.get_close_matches(
        user_disease_n,
        [normalize_text(x) for x in disease_lookup],
        n=5,
        cutoff=0.75,
    )
    if close_matches:
        for original in disease_lookup:
            if normalize_text(original) == close_matches[0]:
                return original

    prompt = f"""
You are a biomedical disease name normalizer.

Task:
Correct spelling mistakes in the user disease term.
Return ONLY valid JSON:
{{
  "corrected_disease": "..."
}}

Examples:
User input: "sytmeic lupus erihtematosus"
Output:
{{"corrected_disease":"systemic lupus erythematosus"}}

User input: "chrnoic lymhocytic leukemia"
Output:
{{"corrected_disease":"chronic lymphocytic leukemia"}}

User input: "{user_disease}"
Output:
""".strip()

    try:
        raw = llm_text(
            prompt, client=client, model_name=model_name,
            max_output_tokens=120,
        )
        parsed = json.loads(raw)
        corrected = parsed.get("corrected_disease", "").strip()
        return corrected if corrected else user_disease
    except Exception:
        if logger:
            logger.warning(
                "Disease name correction failed; using original term."
            )
        return user_disease


def generate_disease_variants(
    user_disease: str,
    client: OpenAI,
    model_name: str = MODEL_NAME,
    logger: Optional[logging.Logger] = None,
) -> List[str]:
    """Generate synonyms / variant spellings of a disease name via LLM.

    Args:
        user_disease: Canonical disease name.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        logger: Optional logger instance.

    Returns:
        Up to six variant names (always including *user_disease* first).
    """
    prompt = f"""
You are a biomedical terminology assistant.

Task:
Generate 5 useful disease name variants or synonyms for database matching.
Keep them highly relevant and not too broad.
Return ONLY valid JSON in this format:
{{
  "variants": ["v1", "v2", "v3", "v4", "v5"]
}}

Examples:
Input disease: lupus
Output:
{{
  "variants": [
    "systemic lupus erythematosus",
    "lupus erythematosus",
    "cutaneous lupus erythematosus",
    "lupus nephritis",
    "systemic autoimmune lupus"
  ]
}}

Input disease: {user_disease}
Output:
""".strip()

    try:
        raw = llm_text(
            prompt, client=client, model_name=model_name,
            max_output_tokens=220,
        )
        parsed = json.loads(raw)
        variants = parsed.get("variants", [])
        variants = [v.strip() for v in variants if str(v).strip()]
        variants = list(dict.fromkeys([user_disease] + variants))
        return variants[:6]
    except Exception:
        if logger:
            logger.warning(
                "Disease variant generation failed; using original term only."
            )
        return [user_disease]


def _score_single_row(
    pos: int,
    row: pd.Series,
    disease_input: str,
    biosample_input: str,
    client: OpenAI,
    model_name: str,
) -> Tuple[int, pd.Series, int, str, Optional[Exception]]:
    """Score one GWAS row against disease + tissue (runs in a worker thread)."""
    trait = str(row.get("DISEASE/TRAIT", ""))
    study = str(row.get("STUDY", ""))

    prompt = f"""
You are a biomedical expert.

Score how well this GWAS row matches:

Disease: {disease_input}
Tissue: {biosample_input}

Scoring rules:
- Disease exact match -> +50
- Disease subtype match -> +40
- Tissue explicitly mentioned -> +40
- Tissue indirectly related -> +20
- No tissue info -> 0
- Wrong tissue -> -30
- Mixed disease (MTAG etc) -> -50

Return ONLY valid JSON:
{{
  "score": 0-100,
  "reason": "short explanation"
}}

Trait: {trait}
Study: {study}
""".strip()

    try:
        raw = llm_text(
            prompt, client=client, model_name=model_name,
            max_output_tokens=150,
        )
        parsed = json.loads(raw)
        score = max(0, min(100, int(parsed.get("score", 0))))
        reason = str(parsed.get("reason", "")).strip()
        return pos, row, score, reason, None
    except Exception as e:
        return pos, row, 0, f"LLM scoring failed: {e}", e


def llm_rank_rows(
    candidate_df: pd.DataFrame,
    disease_input: str,
    biosample_input: str,
    client: OpenAI,
    model_name: str = MODEL_NAME,
    logger: Optional[logging.Logger] = None,
    progress_every: int = STEP6_PROGRESS_EVERY,
    max_workers: int = 10,
) -> pd.DataFrame:
    """Score each GWAS candidate row for disease + tissue relevance via LLM.

    Rows are scored **concurrently** using a thread pool (controlled by
    *max_workers*) for significantly faster execution.

    Two columns are appended: ``LLM_SCORE`` (0--100) and ``LLM_REASON``.

    Args:
        candidate_df: Pre-filtered GWAS rows to rank.
        disease_input: Corrected disease name.
        biosample_input: Normalised tissue / biosample term.
        client: Initialised OpenAI client.
        model_name: Model identifier.
        logger: Optional logger instance.
        progress_every: Log progress every *N* completed rows.
        max_workers: Number of concurrent LLM threads.

    Returns:
        DataFrame sorted by ``LLM_SCORE`` descending.
    """
    total_rows = len(candidate_df)
    if logger:
        logger.info(
            "Ranking %d GWAS rows with LLM scoring (%d concurrent workers)...",
            total_rows, max_workers,
        )

    work_items = [
        (pos, row)
        for pos, (_, row) in enumerate(candidate_df.iterrows(), start=1)
    ]

    ranked_rows = []
    completed = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _score_single_row,
                pos, row, disease_input, biosample_input, client, model_name,
            ): pos
            for pos, row in work_items
        }

        for future in as_completed(futures):
            pos, row, score, reason, error = future.result()
            completed += 1

            row_copy = row.copy()
            row_copy["LLM_SCORE"] = score
            row_copy["LLM_REASON"] = reason
            ranked_rows.append(row_copy)

            if error and logger:
                logger.warning(
                    "Scoring fallback at row %d/%d", pos, total_rows,
                )

            if logger and (
                completed == 1
                or completed % progress_every == 0
                or completed == total_rows
            ):
                logger.info(
                    "  Ranked %d/%d rows", completed, total_rows,
                )

    if not ranked_rows:
        return pd.DataFrame(
            columns=list(candidate_df.columns) + ["LLM_SCORE", "LLM_REASON"]
        )

    ranked_df = pd.DataFrame(ranked_rows)
    return ranked_df.sort_values(
        by="LLM_SCORE", ascending=False
    ).reset_index(drop=True)
