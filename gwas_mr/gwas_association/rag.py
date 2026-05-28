"""GWAS trait retrieval-augmented index: build, cache, and load."""

import hashlib
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .helpers import get_trait_match_column, normalize_text, tokenize_for_index


def get_rag_index_cache_path(index_file: Path, rag_cache_dir: Path) -> Path:
    """Derive a deterministic cache filename from the GWAS index file.

    The hash incorporates the file name, size, and modification time so
    that the cache is automatically invalidated when the source changes.

    Args:
        index_file: Path to the GWAS association TSV.
        rag_cache_dir: Directory where cache files are stored.

    Returns:
        Full path to the JSON cache file.
    """
    stat = index_file.stat()
    cache_key = f"{index_file.name}|{stat.st_size}|{int(stat.st_mtime)}"
    cache_hash = hashlib.md5(cache_key.encode("utf-8")).hexdigest()[:12]
    return rag_cache_dir / f"gwas_trait_rag_index_{cache_hash}.json"


def build_or_load_gwas_rag_index(
    gwas_df: pd.DataFrame,
    index_file: Path,
    rag_cache_dir: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, object]:
    """Build or load a cached GWAS trait RAG index.

    The index maps normalised trait strings to row indices and maintains a
    token-to-trait inverted index for fast lookup.

    Args:
        gwas_df: Loaded GWAS association DataFrame.
        index_file: Path to the source GWAS TSV (used for cache-key hashing).
        rag_cache_dir: Directory for the cache file.  Defaults to
            ``index_file.parent / "_rag_index"``.
        logger: Optional logger instance.

    Returns:
        Dictionary with keys ``trait_column``, ``normalized_traits``,
        ``trait_to_rows``, ``token_to_traits``, and ``unique_traits``.
    """
    if rag_cache_dir is None:
        rag_cache_dir = index_file.parent / "_rag_index"
    rag_cache_dir.mkdir(parents=True, exist_ok=True)

    cache_path = get_rag_index_cache_path(index_file, rag_cache_dir)
    trait_col = get_trait_match_column(gwas_df)

    if cache_path.exists():
        if logger:
            logger.info("Using cached GWAS retrieval index: %s", cache_path.name)
        with open(cache_path, "r", encoding="utf-8") as f:
            rag_index = json.load(f)

        cached_trait_col = rag_index.get("trait_column")
        if cached_trait_col and cached_trait_col in gwas_df.columns:
            return rag_index

        if logger:
            logger.warning(
                "Cached index has invalid trait column — rebuilding."
            )

    if logger:
        logger.info("Building GWAS retrieval index for faster disease lookup...")

    normalized_traits = (
        gwas_df[trait_col]
        .fillna("")
        .astype(str)
        .apply(normalize_text)
        .tolist()
    )

    trait_to_rows: Dict[str, List[int]] = {}
    token_to_traits: Dict[str, List[str]] = {}

    for row_idx, trait in enumerate(normalized_traits):
        if not trait:
            continue
        trait_to_rows.setdefault(trait, []).append(row_idx)

    for trait in trait_to_rows:
        for token in tokenize_for_index(trait):
            token_to_traits.setdefault(token, []).append(trait)

    rag_index = {
        "trait_column": trait_col,
        "normalized_traits": normalized_traits,
        "trait_to_rows": trait_to_rows,
        "token_to_traits": token_to_traits,
        "unique_traits": sorted(trait_to_rows.keys()),
    }

    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(rag_index, f)

    if logger:
        logger.info("Saved GWAS retrieval index: %s", cache_path.name)

    return rag_index


def build_gwas_rag_index_only(
    gwas_index_file: Path,
    rag_cache_dir: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Path]:
    """Convenience wrapper: load the GWAS TSV and build / cache the RAG index.

    Args:
        gwas_index_file: Path to the GWAS association TSV file.
        rag_cache_dir: Directory for the cache file.  Defaults to
            ``gwas_index_file.parent / "_rag_index"``.
        logger: Optional logger instance.

    Returns:
        Dict with ``cache_file`` and ``gwas_index_file`` paths.

    Raises:
        FileNotFoundError: If *gwas_index_file* does not exist.
        ValueError: If the required trait column is missing.
    """
    if not gwas_index_file.exists():
        raise FileNotFoundError(
            f"GWAS index file not found: {gwas_index_file}"
        )

    if logger:
        logger.info("Loading GWAS association file: %s", gwas_index_file)

    gwas_df = pd.read_csv(gwas_index_file, sep="\t", low_memory=False)

    if logger:
        logger.info("Loaded GWAS rows: %s", f"{len(gwas_df):,}")

    trait_col = get_trait_match_column(gwas_df)
    if trait_col not in gwas_df.columns:
        raise ValueError("Missing required trait column in GWAS index")

    if rag_cache_dir is None:
        rag_cache_dir = gwas_index_file.parent / "_rag_index"

    rag_index = build_or_load_gwas_rag_index(
        gwas_df,
        gwas_index_file,
        rag_cache_dir=rag_cache_dir,
        logger=logger,
    )

    cache_path = get_rag_index_cache_path(gwas_index_file, rag_cache_dir)

    if logger:
        logger.info("Index build complete.")
        logger.info(
            "Trait column: %s | Unique traits: %s",
            rag_index.get("trait_column", trait_col),
            f"{len(rag_index.get('unique_traits', [])):,}",
        )

    return {"cache_file": cache_path, "gwas_index_file": gwas_index_file}
