"""gwas_mr -- GWAS & eQTL retrieval and Mendelian Randomization pipeline.

Quick start
-----------
.. code-block:: python

    from gwas_mr import run_full_pipeline

    results = run_full_pipeline(
        disease_name="Breast Cancer",
        biosample_type="Breast Mammary Tissue",
        output_dir="./pipeline_output",
        gwas_data_dir="./gwas_data",
    )

Reference data files (GWAS catalog TSV, GTEx eQTLs, GTEx biosample CSV)
default to the ``gwas_mr_reference/`` directory at the repository root.
The ``mr_pipeline.R`` script is bundled inside the ``gwas_mr`` package.

Individual stages can also be imported separately::

    from gwas_mr import retrieve_gwas_data   # GWAS retrieval only
    from gwas_mr import run_mr               # MR execution only
"""

from .pipeline import run_full_pipeline, preflight_check
from .retrieval import retrieve_gwas_data
from .mr_runner import run_mr
from .report import generate_report

__all__ = [
    "run_full_pipeline",
    "preflight_check",
    "retrieve_gwas_data",
    "run_mr",
    "generate_report",
]
