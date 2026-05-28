"""Default paths for data files and the bundled R script.

The R script lives inside the ``gwas_mr`` package. Large reference
data files (GWAS catalog TSV, GTEx eQTLs, PLINK panel) are looked up
in several locations, with the first valid match winning:

1. ``../shared_reference/gwas_mr_reference/`` for shared sibling repos.
2. ``shared_reference/gwas_mr_reference/`` inside this repository.
3. ``gwas_mr_reference/`` as the legacy wrapper folder.
4. The repository root itself, when the reference data sits beside
   ``main.py``.
"""

from pathlib import Path

# gwas_mr package directory (contains mr_pipeline.R)
PKG_DIR: Path = Path(__file__).resolve().parent

# Repository root (parent of the gwas_mr package)
REPO_ROOT: Path = PKG_DIR.parent


def _looks_like_reference_dir(path: Path) -> bool:
    """Return True when *path* contains the expected reference assets."""
    markers = [
        path / "GWAS-DATABASE-FTP.tsv",
        path / "GTEx_eQTL_TISSUE_EXPRESSION",
        path / "GTEx_EQTL_SAMPLE_COUNTS",
        path / "plink",
    ]
    return any(marker.exists() for marker in markers)


def _resolve_reference_dir() -> Path:
    """Pick the first directory that looks like a valid reference root."""
    candidates = [
        REPO_ROOT.parent / "shared_reference" / "gwas_mr_reference",
        REPO_ROOT / "shared_reference" / "gwas_mr_reference",
        REPO_ROOT / "gwas_mr_reference",
        REPO_ROOT,
    ]
    for candidate in candidates:
        if candidate.is_dir() and _looks_like_reference_dir(candidate):
            return candidate
    return REPO_ROOT / "gwas_mr_reference"


REF_DATA_DIR: Path = _resolve_reference_dir()

DEFAULT_R_SCRIPT: str = str(PKG_DIR / "mr_pipeline.R")

DEFAULT_TSV_PATH: str = str(REF_DATA_DIR / "GWAS-DATABASE-FTP.tsv")

DEFAULT_EQTL_ROOT: str = str(REF_DATA_DIR / "GTEx_eQTL_TISSUE_EXPRESSION")

DEFAULT_GTEX_BIOSAMPLE_CSV: str = str(
    REF_DATA_DIR / "GTEx_EQTL_SAMPLE_COUNTS" / "EQTL-SAMPLE-COUNT-GTEx Portal.csv"
)

DEFAULT_PLINK_DIR: Path = REF_DATA_DIR / "plink"

DEFAULT_GWAS_ASSOC_INDEX: str = str(REF_DATA_DIR / "gwas-latest-all-assoc.tsv")
