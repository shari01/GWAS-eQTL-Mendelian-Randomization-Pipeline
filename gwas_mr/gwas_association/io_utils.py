"""Input file detection and tabular-data loading utilities."""

from pathlib import Path
from typing import Optional

import pandas as pd


def read_table_auto(file_path: Path) -> pd.DataFrame:
    """Read a tabular file, auto-detecting the format from its extension.

    Supported formats: ``.csv``, ``.tsv``, ``.txt`` (tab -> comma ->
    whitespace fallback), ``.xlsx``, and ``.xls``.

    Args:
        file_path: Path to the data file.

    Returns:
        The loaded DataFrame.

    Raises:
        ValueError: If the file extension is not recognised.
    """
    suffix = file_path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(file_path, low_memory=False)
    if suffix == ".tsv":
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    if suffix == ".txt":
        try:
            return pd.read_csv(file_path, sep="\t", low_memory=False)
        except Exception:
            try:
                return pd.read_csv(file_path, sep=",", low_memory=False)
            except Exception:
                return pd.read_csv(file_path, sep=r"\s+", low_memory=False)
    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(file_path)

    raise ValueError(f"Unsupported file type: {file_path.name}")


def detect_input_file(input_folder: Path) -> Path:
    """Auto-detect the best input file in *input_folder*.

    When multiple candidates exist, files whose names contain terms like
    *patient*, *transcript*, *gene*, *deg*, or *expression* are preferred.

    Args:
        input_folder: Directory to scan.

    Returns:
        Path to the chosen input file.

    Raises:
        FileNotFoundError: If no suitable file is found.
    """
    allowed = [".csv", ".tsv", ".txt", ".xlsx", ".xls"]
    files = [
        f for f in input_folder.iterdir()
        if f.is_file() and f.suffix.lower() in allowed
    ]

    if not files:
        raise FileNotFoundError(
            f"No CSV/TSV/TXT/Excel file found in {input_folder}"
        )

    if len(files) == 1:
        return files[0]

    priority_names = ["patient", "transcript", "gene", "deg", "expression"]
    ranked = sorted(
        files,
        key=lambda f: (
            -sum(1 for p in priority_names if p in f.name.lower()),
            f.name.lower(),
        ),
    )
    return ranked[0]


def resolve_input_file(
    input_folder: Path,
    requested_file_name: Optional[str] = None,
) -> Path:
    """Resolve an input file by name, falling back to auto-detection.

    Args:
        input_folder: Directory containing candidate files.
        requested_file_name: Explicit filename to look for (case-insensitive).
            When *None* or blank, :func:`detect_input_file` is used.

    Returns:
        Path to the resolved input file.

    Raises:
        FileNotFoundError: If the requested file does not exist and no
            fallback is available.
    """
    if requested_file_name and str(requested_file_name).strip():
        requested_name = str(requested_file_name).strip()
        requested_path = input_folder / requested_name

        if requested_path.exists() and requested_path.is_file():
            return requested_path

        allowed = [".csv", ".tsv", ".txt", ".xlsx", ".xls"]
        candidates = [
            f for f in input_folder.iterdir()
            if f.is_file() and f.suffix.lower() in allowed
        ]

        for candidate in candidates:
            if candidate.name.lower() == requested_name.lower():
                return candidate

        available_files = ", ".join(sorted(f.name for f in candidates))
        raise FileNotFoundError(
            f"Requested input file '{requested_name}' was not found in "
            f"{input_folder}.  Available files: {available_files}"
        )

    return detect_input_file(input_folder)
