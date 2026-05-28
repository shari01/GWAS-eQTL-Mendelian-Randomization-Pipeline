"""Pipeline configuration constants and lazy OpenAI client factory.

This module provides immutable constants (gene-column candidates, biosample
hints, model defaults) and a lazy-init factory for the OpenAI client so that
**no side effects occur at import time**.
"""

import os
from typing import Dict, List, Optional

# ── Model defaults ──────────────────────────────────────────────────────────

MODEL_NAME: str = "gpt-5.4"
MAX_STUDY_ROWS_FOR_LLM: int = 300
STEP6_PROGRESS_EVERY: int = 10

# ── Gene column detection candidates ────────────────────────────────────────

GENE_COLUMN_CANDIDATES: List[str] = [
    "gene",
    "genes",
    "gene_symbol",
    "gene symbol",
    "gene_symbols",
    "symbol",
    "hgnc_symbol",
    "Gene",
    "Gene_Symbol",
    "Gene Symbol",
]

# ── Biosample / tissue normalisation hints ──────────────────────────────────

TEXTUAL_BIOSAMPLE_HINTS: Dict[str, List[str]] = {
    "pbmc": [
        "pbmc", "peripheral blood mononuclear cell",
        "peripheral blood mononuclear cells",
    ],
    "whole blood": ["whole blood", "blood", "peripheral blood"],
    "blood": ["blood", "whole blood", "peripheral blood", "serum", "plasma"],
    "biopsy": ["biopsy", "tumor biopsy", "tissue biopsy"],
    "tissue": ["tissue", "tumor tissue", "lesional tissue", "organ tissue"],
    "skin": ["skin", "cutaneous", "dermal", "epidermis"],
    "kidney": ["kidney", "renal", "glomerular", "nephritis"],
    "lung": ["lung", "pulmonary", "bronchial", "alveolar"],
    "serum": ["serum", "blood", "plasma"],
    "plasma": ["plasma", "blood", "serum"],
}

# ── Lazy OpenAI client ──────────────────────────────────────────────────────

_openai_client = None


def get_openai_client(api_key: Optional[str] = None) -> "OpenAI":  # noqa: F821
    """Return a shared OpenAI client, creating it on first call.

    The client is cached as a module-level singleton.  Passing an explicit
    *api_key* always creates a fresh client (and replaces the cached one
    only when called without a key next time).

    Args:
        api_key: Explicit API key.  Falls back to the ``OPENAI_API_KEY``
            environment variable (loaded via *python-dotenv* if present).

    Returns:
        An initialised ``openai.OpenAI`` client instance.

    Raises:
        ValueError: If no API key can be resolved.
    """
    global _openai_client

    if _openai_client is not None and api_key is None:
        return _openai_client

    from dotenv import load_dotenv
    from openai import OpenAI

    load_dotenv()

    resolved_key = api_key or os.getenv("OPENAI_API_KEY")
    if not resolved_key:
        raise ValueError(
            "OPENAI_API_KEY not found.  Pass api_key= or set the "
            "OPENAI_API_KEY environment variable."
        )

    client = OpenAI(api_key=resolved_key)

    if api_key is None:
        _openai_client = client

    return client
