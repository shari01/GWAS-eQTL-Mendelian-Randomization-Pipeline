"""Mendelian Randomization pipeline execution via Rscript subprocess.

Wraps the call to ``mr_pipeline.R`` so it can be invoked programmatically
instead of through a CLI.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional


def _log(msg: str) -> None:
    print(msg, flush=True)


def _find_rscript() -> Optional[str]:
    """Locate the ``Rscript`` executable.

    Search order:
    1. System PATH (``shutil.which``)
    2. ``R_HOME`` environment variable — checks ``bin/Rscript``,
       ``bin/x64/Rscript``, and ``bin/Rscript.exe`` (Windows).
    """
    found = shutil.which("Rscript")
    if found:
        return found

    r_home = os.environ.get("R_HOME")
    if not r_home:
        return None

    r_home = Path(r_home)
    candidates = [
        r_home / "bin" / "Rscript.exe",
        r_home / "bin" / "x64" / "Rscript.exe",
        r_home / "bin" / "Rscript",
        r_home / "bin" / "x64" / "Rscript",
    ]
    for c in candidates:
        if c.is_file():
            _log(f"[MR] Found Rscript via R_HOME: {c}")
            return str(c)

    return None


def find_r_script(r_script_path: Optional[str] = None) -> Path:
    """Resolve the location of ``mr_pipeline.R``.

    Search order when *r_script_path* is not provided:
    1. Bundled copy inside ``gwas_mr/`` package
    2. Current working directory
    3. Repository root
    """
    if r_script_path:
        p = Path(r_script_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"R script not found: {r_script_path}")

    from .defaults import DEFAULT_R_SCRIPT

    candidates = [
        Path(DEFAULT_R_SCRIPT),
        Path.cwd() / "mr_pipeline.R",
        Path(__file__).resolve().parent.parent / "mr_pipeline.R",
    ]
    for c in candidates:
        if c.exists():
            return c

    raise FileNotFoundError(
        "mr_pipeline.R not found. It should be inside the gwas_mr/ package, "
        "the working directory, or the repository root — or pass "
        "r_script_path explicitly."
    )


def _find_plink(plink_bin: Optional[str] = None) -> Optional[str]:
    """Locate the ``plink`` executable.

    Search order:
    1. Explicit *plink_bin* argument
    2. ``PLINK_BIN`` environment variable
    3. System PATH
    4. ``gwas_mr_reference/plink/`` directory at the repo root
    """
    if plink_bin and Path(plink_bin).is_file():
        return str(plink_bin)

    from_env = os.environ.get("PLINK_BIN", "")
    if from_env and Path(from_env).is_file():
        return from_env

    from_path = shutil.which("plink")
    if from_path:
        return from_path

    from .defaults import DEFAULT_PLINK_DIR
    for name in ("plink.exe", "plink"):
        candidate = DEFAULT_PLINK_DIR / name
        if candidate.is_file():
            return str(candidate)

    return None


def _resolve_plink_ref(plink_ref: Optional[str] = None) -> str:
    """Resolve the PLINK reference panel prefix.

    Search order:
    1. Explicit *plink_ref* argument
    2. ``PLINK_REF`` environment variable
    3. ``gwas_mr_reference/plink/1000G_EUR_hg38``
    """
    if plink_ref:
        return plink_ref

    from_env = os.environ.get("PLINK_REF", "")
    if from_env:
        return from_env

    from .defaults import DEFAULT_PLINK_DIR
    return str(DEFAULT_PLINK_DIR / "1000G_EUR_hg38")


def run_mr(
    base_dir: str,
    eqtl_path: str,
    gwas_path: str,
    n_eqtl: int,
    n_gwas: int,
    out_dir: str,
    disease_name: str,
    biosample_type: str,
    *,
    r_script_path: Optional[str] = None,
    plink_bin: Optional[str] = None,
    plink_ref: Optional[str] = None,
) -> None:
    """
    Run the Mendelian Randomization pipeline via ``Rscript``.

    Parameters
    ----------
    base_dir : str
        Base directory used by the R pipeline for config context.
    eqtl_path : str
        Absolute path to the eQTL file.
    gwas_path : str
        Absolute path to the GWAS harmonised TSV file.
    n_eqtl : int
        eQTL sample size.
    n_gwas : int
        GWAS sample size.
    out_dir : str
        Output directory for MR results.
    disease_name : str
        Disease / outcome label used in evidence-table filenames.
    biosample_type : str
        Biosample type label used in evidence-table filenames.
    r_script_path : str, optional
        Path to ``mr_pipeline.R``. Auto-detected if not given.
    plink_bin : str, optional
        Path to the PLINK 1.9 executable.  Auto-detected from
        ``PLINK_BIN`` env var, system PATH, or ``gwas_mr_reference/plink/``.
    plink_ref : str, optional
        Prefix for the PLINK reference panel (``.bed/.bim/.fam``).
        Auto-detected from ``PLINK_REF`` env var or
        ``gwas_mr_reference/plink/1000G_EUR_hg38``.

    Raises
    ------
    FileNotFoundError
        If required paths (base_dir, eqtl_path, gwas_path, R script) are
        missing.
    RuntimeError
        If the Rscript process exits with a non-zero return code.
    """
    rscript_bin = _find_rscript()
    if rscript_bin is None:
        raise RuntimeError(
            "Rscript not found. Checked system PATH and R_HOME. "
            "Either add R's bin/ directory to your PATH, or set R_HOME "
            "to your R installation root (e.g. C:\\Program Files\\R\\R-4.x.x)."
        )

    r_script = find_r_script(r_script_path)

    for label, path in [
        ("base_dir", base_dir),
        ("eqtl_path", eqtl_path),
        ("gwas_path", gwas_path),
    ]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} does not exist: {path}")

    os.makedirs(out_dir, exist_ok=True)

    # Resolve PLINK paths and pass to R script via environment variables
    resolved_plink = _find_plink(plink_bin)
    resolved_ref = _resolve_plink_ref(plink_ref)

    env = os.environ.copy()
    if resolved_plink:
        env["PLINK_BIN"] = resolved_plink
        _log(f"[MR] PLINK binary: {resolved_plink}")
    else:
        _log("[MR] PLINK not found — LD clumping / LD matrix will be skipped.")

    if resolved_ref:
        env["PLINK_REF"] = resolved_ref
        _log(f"[MR] PLINK ref panel: {resolved_ref}")

    cmd = [
        rscript_bin,
        str(r_script),
        str(base_dir),
        str(eqtl_path),
        str(gwas_path),
        str(n_eqtl),
        str(n_gwas),
        str(out_dir),
        str(disease_name),
        str(biosample_type),
    ]

    _log(f"[MR] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, env=env, check=False)

    if result.returncode != 0:
        raise RuntimeError(
            f"MR pipeline failed with exit code {result.returncode}"
        )

    _log(f"[MR] Pipeline completed successfully. Results in: {out_dir}")
