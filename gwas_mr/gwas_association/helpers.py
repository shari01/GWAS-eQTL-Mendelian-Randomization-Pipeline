"""Shared text normalisation and filesystem-safe naming utilities."""

import re
from typing import List

import pandas as pd


def get_trait_match_column(gwas_df: pd.DataFrame) -> str:
    """Return the best available trait column name in a GWAS DataFrame.

    Args:
        gwas_df: GWAS association DataFrame.

    Returns:
        ``'MAPPED_TRAIT'`` if present, otherwise ``'DISEASE/TRAIT'``.
    """
    if "MAPPED_TRAIT" in gwas_df.columns:
        return "MAPPED_TRAIT"
    return "DISEASE/TRAIT"


def normalize_text(x: str) -> str:
    """Lowercase, strip, and collapse internal whitespace.

    Args:
        x: Raw text value (may be NaN-like).

    Returns:
        Cleaned lowercase string, or ``""`` for NaN / None values.
    """
    if pd.isna(x):
        return ""
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x


def safe_filename(text: str) -> str:
    """Sanitise *text* into a safe filesystem name (max 120 chars).

    Args:
        text: Arbitrary string.

    Returns:
        A string containing only word characters, hyphens, and dots.
    """
    text = re.sub(r"[^\w\-\.]+", "_", text.strip())
    text = re.sub(r"_+", "_", text)
    return text[:120]


def clean_name(text: str) -> str:
    """Return a lowercase, filesystem-safe version of *text*.

    Args:
        text: Arbitrary string.

    Returns:
        Sanitised lowercase string.
    """
    return safe_filename(str(text).lower())


def tokenize_for_index(text: str) -> List[str]:
    """Split normalised text into alphanumeric tokens of length >= 2.

    Args:
        text: Raw text to tokenise.

    Returns:
        List of lowercase tokens suitable for inverted-index lookup.
    """
    normalized = normalize_text(text)
    tokens = re.findall(r"[a-z0-9]+", normalized)
    return [tok for tok in tokens if len(tok) >= 2]
