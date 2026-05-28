"""End-to-end GWAS retrieval + Mendelian Randomization orchestrator.

The single entry-point :func:`run_full_pipeline` chains GWAS data
retrieval with MR analysis and writes all results under one output
directory.
"""

import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from .retrieval import log, retrieve_gwas_data, safe_name
from .mr_runner import run_mr, _find_rscript, find_r_script, _find_plink, _resolve_plink_ref
from .report import generate_report


# ── Post-run output cleanup (safety net) ───────────────────────────────────
#
# The R pipeline already writes per-gene intermediates to tempdir and
# skips plots for <2-SNP genes.  This function acts as a lightweight
# fallback that removes any leftover empty directories.

_EPHEMERAL_DIRS = [
    Path("03_instruments", "03a_eqtl_filtered", "per_gene"),
    Path("05_clumping"),
    Path("06_harmonisation"),
    Path("07_ld"),
    Path("08_mr"),
    Path("09_finemap"),
    Path("10_coloc"),
    Path("11_plots", "summary"),
]


def _cleanup_output(mr_output_dir: str) -> None:
    """Remove any leftover empty directories from the MR output."""
    out = Path(mr_output_dir)
    if not out.is_dir():
        return

    for rel in _EPHEMERAL_DIRS:
        d = out / rel
        if d.is_dir():
            try:
                shutil.rmtree(d)
            except OSError:
                pass
            parent = d.parent
            while parent != out and parent.is_dir() and not any(parent.iterdir()):
                parent.rmdir()
                parent = parent.parent

    plots_dir = out / "11_plots" / "per_gene"
    if plots_dir.is_dir() and not any(plots_dir.iterdir()):
        plots_dir.rmdir()
    plots_parent = out / "11_plots"
    if plots_parent.is_dir() and not any(plots_parent.iterdir()):
        plots_parent.rmdir()


# ── Preflight dependency check ─────────────────────────────────────────────

def _check_python_imports() -> List[str]:
    """Return a list of missing Python packages required by the pipeline."""
    missing = []
    required = {
        "pandas": "pandas",
        "requests": "requests",
        "openpyxl": "openpyxl",
    }
    for import_name, pip_name in required.items():
        try:
            __import__(import_name)
        except ImportError:
            missing.append(pip_name)
    return missing


def _check_openai_key() -> bool:
    """Return True if the OpenAI API key is discoverable."""
    try:
        from decouple import config, UndefinedValueError
        config("OPENAI_API_KEY")
        return True
    except Exception:
        return False


def preflight_check(
    *,
    tsv_path: Optional[str] = None,
    eqtl_root: Optional[str] = None,
    gtex_biosample_csv: Optional[str] = None,
    r_script_path: Optional[str] = None,
    plink_bin: Optional[str] = None,
    plink_ref: Optional[str] = None,
    use_llm: bool = True,
    run_mr_analysis: bool = True,
) -> bool:
    """Validate that all dependencies are in place before running.

    Prints a status report and returns ``True`` if everything required
    is present.  Missing *optional* items produce warnings but do not
    cause the check to fail.
    """
    from .defaults import (
        DEFAULT_TSV_PATH,
        DEFAULT_EQTL_ROOT,
        DEFAULT_GTEX_BIOSAMPLE_CSV,
        DEFAULT_PLINK_DIR,
        REF_DATA_DIR,
    )

    ok_sym = "[OK]   "
    warn_sym = "[WARN] "
    fail_sym = "[FAIL] "

    all_ok = True

    log("=" * 60)
    log("PREFLIGHT DEPENDENCY CHECK")
    log("=" * 60)

    # ── 1. Python packages ─────────────────────────────────────────────
    log("\n--- Python packages ---")
    missing_py = _check_python_imports()
    if missing_py:
        all_ok = False
        log(f"{fail_sym}Missing Python packages: {', '.join(missing_py)}")
        log(f"       Fix: pip install {' '.join(missing_py)}")
    else:
        log(f"{ok_sym}All required Python packages installed (pandas, requests, openpyxl)")

    if use_llm:
        try:
            __import__("openai")
            log(f"{ok_sym}openai package installed")
        except ImportError:
            all_ok = False
            log(f"{fail_sym}openai package not installed (needed when use_llm=True)")
            log("       Fix: pip install openai")

        try:
            __import__("decouple")
            log(f"{ok_sym}python-decouple package installed")
        except ImportError:
            all_ok = False
            log(f"{fail_sym}python-decouple package not installed")
            log("       Fix: pip install python-decouple")

    # ── 2. OpenAI API key ──────────────────────────────────────────────
    log("\n--- OpenAI API key ---")
    if use_llm:
        if _check_openai_key():
            log(f"{ok_sym}OPENAI_API_KEY found")
        else:
            all_ok = False
            log(f"{fail_sym}OPENAI_API_KEY not found")
            log("       Fix: Create a .env file with OPENAI_API_KEY=sk-...")
            log("            or set it as a system environment variable.")
    else:
        log(f"{ok_sym}LLM disabled (use_llm=False) — API key not required")

    # ── 3. Reference data files ────────────────────────────────────────
    log("\n--- Reference data files ---")
    _tsv = tsv_path or DEFAULT_TSV_PATH
    _eqtl = eqtl_root or DEFAULT_EQTL_ROOT
    _gtex = gtex_biosample_csv or DEFAULT_GTEX_BIOSAMPLE_CSV

    ref_label = str(REF_DATA_DIR)

    log(f"  Reference dir: {ref_label}")

    if Path(_tsv).is_file():
        log(f"{ok_sym}GWAS catalog TSV: {_tsv}")
    else:
        all_ok = False
        log(f"{fail_sym}GWAS catalog TSV not found: {_tsv}")
        log(f"       Fix: Download from the Google Drive link and place in {ref_label}/")

    if Path(_eqtl).is_dir():
        biosample_count = len(list(Path(_eqtl).glob("*.txt")))
        log(f"{ok_sym}GTEx eQTL directory: {_eqtl} ({biosample_count} biosample files)")
    else:
        all_ok = False
        log(f"{fail_sym}GTEx eQTL directory not found: {_eqtl}")
        log(f"       Fix: Download from the Google Drive link and extract into {ref_label}/")

    if Path(_gtex).is_file():
        log(f"{ok_sym}GTEx biosample CSV: {_gtex}")
    else:
        all_ok = False
        log(f"{fail_sym}GTEx biosample sample-count CSV not found: {_gtex}")
        log(f"       Fix: Download from the Google Drive link and extract into {ref_label}/")

    # ── 4. R / Rscript (only if MR is enabled) ─────────────────────────
    if run_mr_analysis:
        log("\n--- R / Rscript ---")
        rscript = _find_rscript()
        if rscript:
            log(f"{ok_sym}Rscript found: {rscript}")
        else:
            all_ok = False
            log(f"{fail_sym}Rscript not found on PATH or via R_HOME")
            log("       Fix: Install R and add to PATH, or set R_HOME env var")

        # ── 5. mr_pipeline.R ───────────────────────────────────────────
        log("\n--- mr_pipeline.R ---")
        try:
            r_script = find_r_script(r_script_path)
            log(f"{ok_sym}R script: {r_script}")
        except FileNotFoundError:
            all_ok = False
            log(f"{fail_sym}mr_pipeline.R not found")
            log("       Fix: Ensure it exists in the gwas_mr/ package directory")

        # ── 6. PLINK ───────────────────────────────────────────────────
        log("\n--- PLINK 1.9 ---")
        resolved_plink = _find_plink(plink_bin)
        if resolved_plink:
            log(f"{ok_sym}PLINK binary: {resolved_plink}")
        else:
            all_ok = False
            log(f"{fail_sym}PLINK not found")
            log("       Fix: Download from https://www.cog-genomics.org/plink/")
            log(f"            Place in {ref_label}/plink/ or set PLINK_BIN env var")

        resolved_ref = _resolve_plink_ref(plink_ref)
        ref_bed = Path(resolved_ref + ".bed")
        if ref_bed.is_file():
            log(f"{ok_sym}PLINK reference panel: {resolved_ref}")
        else:
            all_ok = False
            log(f"{fail_sym}PLINK reference panel not found: {resolved_ref}")
            log("       Fix: Download 1000G EUR panel and place .bed/.bim/.fam in")
            log(f"            {ref_label}/plink/ as 1000G_EUR_hg38.*")

        # ── 7. R packages (check via Rscript) ─────────────────────────
        if rscript:
            log("\n--- R packages ---")
            r_check_code = (
                'all_pkgs <- c("data.table","TwoSampleMR","ggplot2",'
                '"coloc","susieR","pheatmap","jsonlite");'
                'missing <- all_pkgs[!sapply(all_pkgs,requireNamespace,quietly=TRUE)];'
                'if(length(missing)) cat("MISSING:",paste(missing,collapse=","),"\\n") '
                'else cat("ALL_OK\\n")'
            )
            import subprocess
            try:
                proc = subprocess.run(
                    [rscript, "-e", r_check_code],
                    capture_output=True, text=True, timeout=60,
                )
                r_out = proc.stdout.strip()
                for line in r_out.splitlines():
                    if line.startswith("MISSING:"):
                        pkgs = line.split(":", 1)[1].strip()
                        all_ok = False
                        log(f"{fail_sym}Missing R packages: {pkgs}")
                        log("       Fix: In R console run:")
                        log('         install.packages(c("data.table","ggplot2","coloc",')
                        log('                           "susieR","pheatmap","jsonlite"))')
                        log('         install.packages("remotes")')
                        log('         remotes::install_github("MRCIEU/TwoSampleMR")')
                    elif line == "ALL_OK":
                        log(f"{ok_sym}All R packages installed (data.table, TwoSampleMR,")
                        log("         ggplot2, coloc, susieR, pheatmap, jsonlite)")
            except (subprocess.TimeoutExpired, FileNotFoundError):
                all_ok = False
                log(f"{fail_sym}Could not verify R packages (Rscript timed out or failed)")
    else:
        log("\n--- R / PLINK ---")
        log(f"{ok_sym}MR analysis disabled (run_mr_analysis=False) — R/PLINK not required")

    # ── Summary ────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    if all_ok:
        log("PREFLIGHT CHECK PASSED — all required dependencies are in place.")
    else:
        log("PREFLIGHT CHECK FAILED — some required dependencies are missing.")
        log("See [FAIL] items above for details.")
    log("=" * 60 + "\n")

    return all_ok


# ── Main pipeline ───────────────────────────────────────────────────────────

def run_full_pipeline(
    disease_name: str,
    biosample_type: str,
    output_dir: str,
    gwas_data_dir: str,
    *,
    tsv_path: Optional[str] = None,
    eqtl_root: Optional[str] = None,
    gtex_biosample_csv: Optional[str] = None,
    r_script_path: Optional[str] = None,
    plink_bin: Optional[str] = None,
    plink_ref: Optional[str] = None,
    use_llm: bool = True,
    top_k: int = 1,
    run_mr_analysis: bool = True,
    run_association: bool = True,
    input_genes: Optional[Union[List[str], Sequence[str], Path]] = None,
    gwas_assoc_index: Optional[str] = None,
    skip_preflight: bool = False,
) -> List[Dict[str, Any]]:
    """
    End-to-end pipeline: retrieve GWAS data, then run MR analysis.

    Before execution a **preflight check** validates that all required
    dependencies (Python packages, API keys, reference data, R, PLINK)
    are available.  If any required dependency is missing the pipeline
    aborts early with a clear report.  Pass ``skip_preflight=True`` to
    bypass this check.

    All data-path parameters are optional and default to the
    ``gwas_mr_reference/`` directory at the repository root.

    Parameters
    ----------
    disease_name : str
        Disease or phenotype name (e.g. ``"Breast Cancer"``).
    biosample_type : str
        GTEx biosample type (e.g. ``"Whole Blood"``, ``"Pancreas"``, ``"PBMC"``).
    output_dir : str
        Root directory for MR results.  MR output goes to
        ``output_dir/mr_results/<trait>/``.
    gwas_data_dir : str
        Directory for downloaded GWAS summary statistics and cached
        files.  Kept separate from ``output_dir`` so that large GWAS
        downloads can be reused across multiple pipeline runs.
    tsv_path : str, optional
        Path to the GWAS catalog TSV.  Defaults to
        ``gwas_mr_reference/GWAS-DATABASE-FTP.tsv``.
    eqtl_root : str, optional
        Path to the GTEx eQTL biosample expression directory.
        Defaults to ``gwas_mr_reference/GTEx_eQTL_TISSUE_EXPRESSION/``.
    gtex_biosample_csv : str, optional
        Path to the GTEx biosample sample-count CSV.  Defaults to
        ``gwas_mr_reference/GTEx_EQTL_SAMPLE_COUNTS/EQTL-SAMPLE-COUNT-GTEx Portal.csv``.
    r_script_path : str, optional
        Path to ``mr_pipeline.R``.  Defaults to the copy bundled
        inside the ``gwas_mr`` package.
    plink_bin : str, optional
        Path to the PLINK 1.9 executable.  Auto-detected from ``PLINK_BIN``
        env var, system PATH, or ``gwas_mr_reference/plink/``.
    plink_ref : str, optional
        Prefix for the PLINK LD reference panel (``.bed/.bim/.fam``).
        Auto-detected from ``PLINK_REF`` env var or
        ``gwas_mr_reference/plink/1000G_EUR_hg38``.
    use_llm : bool
        Use OpenAI LLM for disease-name normalisation (default ``True``).
        The API key is read via ``python-decouple`` (checks ``.env``,
        environment variables, and ``settings.ini`` automatically).
    top_k : int
        Number of top GWAS rows per trait by sample size (default ``1``).
    run_mr_analysis : bool
        Whether to execute the MR step after retrieval (default ``True``).
    run_association : bool
        Whether to execute the GWAS association mapping step in parallel
        with MR (default ``True``).  Requires ``gwas-latest-all-assoc.tsv``
        in ``gwas_mr_reference/`` and a non-empty ``input_genes`` list.
    input_genes : list[str] or Path, optional
        Gene symbols for the GWAS association sub-pipeline.  Either a
        **list of gene symbols** (e.g. ``["FLT3", "NPM1", "RUNX1"]``)
        or a **Path** to a tabular file whose gene column is
        auto-detected.  Required when ``run_association=True``; ignored
        otherwise.
    gwas_assoc_index : str, optional
        Path to the GWAS association index TSV used by the association
        sub-pipeline.  Defaults to
        ``gwas_mr_reference/gwas-latest-all-assoc.tsv``.
    skip_preflight : bool
        Skip the preflight dependency check (default ``False``).

    Returns
    -------
    list[dict]
        One dict per processed trait.  Each dict contains all retrieval
        metadata plus ``mr_output_dir`` and ``mr_status`` when MR was run.
    """
    # ── Preflight check ────────────────────────────────────────────────
    if not skip_preflight:
        passed = preflight_check(
            tsv_path=tsv_path,
            eqtl_root=eqtl_root,
            gtex_biosample_csv=gtex_biosample_csv,
            r_script_path=r_script_path,
            plink_bin=plink_bin,
            plink_ref=plink_ref,
            use_llm=use_llm,
            run_mr_analysis=run_mr_analysis,
        )
        if not passed:
            raise RuntimeError(
                "Preflight check failed — required dependencies are missing. "
                "See the report above for details and fixes. "
                "Pass skip_preflight=True to bypass this check."
            )

    output_path = Path(output_dir)
    gwas_data_path = Path(gwas_data_dir)

    output_path.mkdir(parents=True, exist_ok=True)
    gwas_data_path.mkdir(parents=True, exist_ok=True)
    mr_results_dir = output_path 

    log("=" * 60)
    log("GWAS-eQTL Full Pipeline")
    log(f"  Disease:    {disease_name}")
    log(f"  Biosample:  {biosample_type}")
    log(f"  GWAS data:  {gwas_data_path}")
    log(f"  MR output:  {output_dir}")
    log("=" * 60)

    # ── Step 1: Retrieve GWAS data ─────────────────────────────────────────
    log("\n[STEP 1/2] Retrieving GWAS data...")

    datasets = retrieve_gwas_data(
        disease_name=disease_name,
        biosample_type=biosample_type,
        tsv_path=tsv_path,
        eqtl_root=eqtl_root,
        gtex_biosample_csv=gtex_biosample_csv,
        out_root=str(gwas_data_path),
        use_llm=use_llm,
        top_k=top_k,
    )

    if not datasets:
        log("\n[DONE] No GWAS datasets retrieved. Pipeline stopped.")
        return datasets

    log(f"\n[OK] Retrieved {len(datasets)} dataset(s).")

    # ── Validate association inputs ─────────────────────────────────────
    has_assoc_genes = bool(run_association and input_genes)
    if run_association and not input_genes:
        log("[WARN] run_association=True but no input_genes provided — "
            "association step will be skipped.")

    # ── Step 2: Run MR + GWAS association in parallel ────────────────────
    both_skipped = not run_mr_analysis and not has_assoc_genes
    if both_skipped:
        log("\n[SKIP] Both MR analysis and association analysis are disabled.")
        for ds in datasets:
            ds["mr_status"] = "skipped"
            ds["mr_output_dir"] = None
        return datasets

    step_parts = []
    if run_mr_analysis:
        step_parts.append("Mendelian Randomization")
    if has_assoc_genes:
        step_parts.append("GWAS Association")
    log(f"\n[STEP 2/2] Running {' + '.join(step_parts)} "
        f"{'in parallel ' if len(step_parts) > 1 else ''}...")

    def _mr_worker() -> None:
        """Execute MR analysis for each retrieved dataset."""
        for ds in datasets:
            trait = ds.get("trait", "unknown")
            gwas_path = (
                ds.get("gwas_path")
                or ds.get("gwas_path_tsv")
                or ds.get("gwas_path_gz")
            )
            eqtl_path = ds.get("eqtl_path")
            n_eqtl = ds.get("n_eqtl_for_mr")
            n_gwas = ds.get("N_gwas")

            if not gwas_path or not eqtl_path:
                log(f"\n[SKIP] MR for '{trait}': missing GWAS or eQTL path.")
                ds["mr_status"] = "skipped_missing_paths"
                ds["mr_output_dir"] = None
                continue

            if not n_eqtl or not n_gwas:
                log(f"\n[SKIP] MR for '{trait}': missing sample size(s).")
                ds["mr_status"] = "skipped_missing_sample_sizes"
                ds["mr_output_dir"] = None
                continue

            mr_out = str(mr_results_dir / safe_name(trait))
            ds["mr_output_dir"] = mr_out

            log(f"\n{'─' * 50}")
            log(f"Running MR for: {trait}")
            log(f"  GWAS path:  {gwas_path}")
            log(f"  eQTL path:  {eqtl_path}")
            log(f"  N_eQTL:     {n_eqtl}")
            log(f"  N_GWAS:     {n_gwas}")
            log(f"  Output:     {mr_out}")

            try:
                run_mr(
                    base_dir=str(output_path),
                    eqtl_path=eqtl_path,
                    gwas_path=gwas_path,
                    n_eqtl=int(n_eqtl),
                    n_gwas=int(n_gwas),
                    out_dir=mr_out,
                    disease_name=disease_name,
                    biosample_type=biosample_type,
                    r_script_path=r_script_path,
                    plink_bin=plink_bin,
                    plink_ref=plink_ref,
                )
                ds["mr_status"] = "success"
                log(f"[OK] MR completed for '{trait}'")
            except Exception as e:
                ds["mr_status"] = f"failed: {e}"
                log(f"[ERROR] MR failed for '{trait}': {e}")

            if ds["mr_status"] == "success":
                try:
                    rpt = generate_report(
                        mr_output_dir=mr_out,
                        disease_name=disease_name,
                        biosample_type=biosample_type,
                        gwas_accession=ds.get("accession"),
                        n_gwas=n_gwas,
                        n_eqtl=n_eqtl,
                    )
                    if rpt:
                        log(f"[OK] Summary report: {rpt}")
                except Exception as rpt_err:
                    log(f"[WARN] Report generation failed: {rpt_err}")

                try:
                    _cleanup_output(mr_out)
                except Exception as cleanup_err:
                    log(f"[WARN] Output cleanup failed: {cleanup_err}")

    def _association_worker():
        """Execute the GWAS association mapping pipeline."""
        from .gwas_association import run_gwas_association
        from .defaults import DEFAULT_GWAS_ASSOC_INDEX

        index_path = Path(gwas_assoc_index or DEFAULT_GWAS_ASSOC_INDEX)
        if not index_path.is_file():
            log(f"[WARN] GWAS association index not found: {index_path} — skipping association step.")
            return None

        assoc_output_dir = output_path / safe_name(disease_name)
        genes_desc = (
            str(input_genes) if isinstance(input_genes, Path)
            else f"{len(input_genes)} gene(s)"
        )
        log(f"\n{'─' * 50}")
        log(f"Running GWAS Association mapping")
        log(f"  Index file:  {index_path}")
        log(f"  Input genes: {genes_desc}")
        log(f"  Output:      {assoc_output_dir / 'gwas_association'}")

        try:
            result = run_gwas_association(
                gwas_index_file=index_path,
                input_genes=input_genes,
                output_dir=assoc_output_dir,
                disease=disease_name,
                biosample=biosample_type,
            )
            log(f"[OK] GWAS Association completed — "
                f"{result.gwas_rows_selected} rows selected, "
                f"{result.genes_mapped} genes mapped.")
            return result
        except Exception as e:
            log(f"[ERROR] GWAS Association failed: {e}")
            return None

    # ── Launch workers ───────────────────────────────────────────────────
    assoc_result = None

    if not run_mr_analysis and has_assoc_genes:
        log("\n[SKIP] MR analysis skipped (run_mr_analysis=False).")
        for ds in datasets:
            ds["mr_status"] = "skipped"
            ds["mr_output_dir"] = None
        assoc_result = _association_worker()

    elif run_mr_analysis and not has_assoc_genes:
        if not run_association:
            log("\n[SKIP] GWAS Association skipped (run_association=False).")
        elif not input_genes:
            log("\n[SKIP] GWAS Association skipped (no input_genes provided).")
        _mr_worker()

    elif run_mr_analysis and has_assoc_genes:
        with ThreadPoolExecutor(max_workers=2) as pool:
            mr_future = pool.submit(_mr_worker)
            assoc_future = pool.submit(_association_worker)

            mr_future.result()
            assoc_result = assoc_future.result()

    # ── Summary ────────────────────────────────────────────────────────────
    succeeded = sum(1 for d in datasets if d.get("mr_status") == "success")
    total = len(datasets)

    log(f"\n{'=' * 60}")
    log(f"Pipeline complete: {succeeded}/{total} trait(s) MR processed successfully.")
    if assoc_result is not None:
        log(f"GWAS Association: {assoc_result.gwas_rows_selected} rows, "
            f"{assoc_result.genes_mapped} genes  →  {assoc_result.run_dir}")
    elif run_association:
        log("GWAS Association: not completed (see warnings above).")
    log(f"Results directory: {output_dir}")
    log("=" * 60)

    return datasets
