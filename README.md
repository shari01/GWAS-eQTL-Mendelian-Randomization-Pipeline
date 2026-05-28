# GWAS Retrieval & Mendelian Randomization Pipeline

An end-to-end Python package that retrieves GWAS summary statistics from the
NHGRI-EBI GWAS Catalog, pairs them with GTEx eQTL data, and runs a
comprehensive Mendelian Randomization (MR) analysis via an R backend.

---

## Repository layout

```
gwas_mrrieval_pipeline/
├── gwas_mr/                       # Python package
│   ├── __init__.py
│   ├── defaults.py                 # Default path resolution
│   ├── retrieval.py                # GWAS data retrieval
│   ├── mr_runner.py                # R subprocess launcher
│   ├── pipeline.py                 # Orchestrator (run_full_pipeline)
│   └── mr_pipeline.R               # Bundled MR R script
├── gwas_mr_reference/              # Large reference data (see step 4)
│   ├── GWAS-DATABASE-FTP.tsv
│   ├── GTEx_eQTL_TISSUE_EXPRESSION/
│   ├── GTEx_EQTL_SAMPLE_COUNTS/
│   └── plink/
├── main.py                         # Example entry point
└── README.md
```

---

## Prerequisites

- **Python 3.9+**
- **R 4.x** (with `Rscript` accessible — see [R setup](#2-install-r-packages))
- **PLINK 1.9** (included in the data download — see [step 3](#3-download-reference-data--plink))

---

## 1. Install Python dependencies

Create a virtual environment and install the required packages:

```bash
python -m venv .venv

# Activate (Windows PowerShell)
.venv\Scripts\Activate.ps1

# Activate (Linux / macOS)
source .venv/bin/activate
```

Then install:

```bash
pip install pandas requests openai python-decouple openpyxl
```

| Package | Purpose |
|---|---|
| `pandas` | Data manipulation (GWAS table filtering, reports) |
| `requests` | HTTP downloads from GWAS Catalog FTP |
| `openai` | LLM-based disease name normalization (optional, see `use_llm` flag) |
| `python-decouple` | Reads `OPENAI_API_KEY` from `.env` / environment variables |
| `openpyxl` | Writing Excel summary reports |

### OpenAI API key (only needed if `use_llm=True`, the default)

Create a `.env` file in your working directory:

```
OPENAI_API_KEY=sk-...your-key-here...
```

`python-decouple` will find it automatically. Alternatively, set it as a system
environment variable.  If you don't have an API key, pass `use_llm=False` to
skip LLM normalization and use your disease name as-is.

---

## 2. Install R packages

Open an R console and run:

```r
# Mandatory (pipeline will not run without these)
install.packages(c("data.table", "ggplot2"))
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("coloc", "susieR", "pheatmap", "jsonlite"))
```

| R Package | What it does |
|---|---|
| `data.table` | All data loading and manipulation |
| `TwoSampleMR` | Core MR methods (IVW, Egger, median, harmonization) |
| `ggplot2` | Scatter, forest, funnel, leave-one-out plots |
| `coloc` | Colocalization analysis (shared causal variant probability) |
| `susieR` | Fine-mapping via SuSiE (credible set identification) |
| `pheatmap` | LD heatmap plots |
| `jsonlite` | Saves run configuration as JSON |

All packages above are **required** — the preflight check will fail if any
are missing.

Make sure `Rscript` is in your system PATH, or set the `R_HOME` environment
variable to your R installation directory (e.g. `C:\Program Files\R\R-4.4.2`).

---

## 3. Download reference data & PLINK

The pipeline requires reference data files and PLINK 1.9 (for LD clumping
and LD matrix computation). Everything is bundled in a single zip on
Google Drive:

**Download link:**
https://drive.google.com/file/d/1MyDiKtGCcsHn-0ZH4pkoup3jlRvba9LM/view?usp=drive_link

Extract the zip contents into the `gwas_mr_reference/` directory at the
repository root. The expected structure is:

```
gwas_mr_reference/
├── GWAS-DATABASE-FTP.tsv                    (GWAS Catalog index with FTP URLs)
├── GTEx_eQTL_TISSUE_EXPRESSION/             (one .txt per biosample type)
│   ├── Breast_Mammary_Tissue.v10.eGenes.txt
│   ├── Whole_Blood.v10.eGenes.txt
│   ├── Brain_Cortex.v10.eGenes.txt
│   └── ... (49 biosample types total)
├── GTEx_EQTL_SAMPLE_COUNTS/
│   └── EQTL-SAMPLE-COUNT-GTEx Portal.csv
└── plink/
    ├── plink.exe                            (or plink on Linux/macOS)
    ├── 1000G_EUR_hg38.bed                   (1000 Genomes EUR LD reference)
    ├── 1000G_EUR_hg38.bim
    └── 1000G_EUR_hg38.fam
```

All files including PLINK and the LD reference panel are included in the
zip. The pipeline auto-detects them from `gwas_mr_reference/plink/`.
Alternatively, set `PLINK_BIN` and `PLINK_REF` environment variables or
pass paths explicitly (see [Advanced configuration](#advanced-configuration)).

---

## 4. Run the pipeline

Create a Python script (e.g. `main.py`) in the repository root:

```python
from gwas_mr import run_full_pipeline

results = run_full_pipeline(
    disease_name="Breast Cancer",
    biosample_type="Breast Mammary Tissue",
    output_dir="pipeline_output",
    gwas_data_dir="gwas_data",
)
```

Then run:

```bash
python main.py
```

### What happens

1. **GWAS retrieval** — The disease name is normalized (via LLM if enabled),
   matched against the GWAS Catalog index, and the best harmonized summary
   statistics are downloaded into `gwas_data/`.
2. **MR analysis** — The R script runs IVW, Wald ratio, MR-Egger, weighted
   median/mode, sensitivity tests, colocalization, and fine-mapping. All
   results, plots, and tables are written to `pipeline_output/`.

### Output structure

```
pipeline_output/<trait>/
├── 00_config/           Run configuration, session info
├── 01_input_raw/        Raw eQTL + GWAS copies
├── 02_input_parsed/     Standardized columns
├── 03_instruments/      Filtered eQTL instruments (F-stat)
├── 04_overlap/          SNPs shared between eQTL and GWAS
├── 05_clumping/         LD clumping results (PLINK)
├── 06_harmonised/       Allele-harmonized exposure/outcome
├── 07_ld/               LD matrix, correlated SNP tables
├── 08_mr/               Core MR results (all methods)
├── 09_diagnostics/      Heterogeneity, pleiotropy, leave-one-out
├── 10_plots/            Scatter, forest, funnel, LD heatmaps
├── 11_coloc_finemap/    Colocalization + SuSiE fine-mapping
└── 12_summary_tables/   Final CSVs: results, QC, interpretation tiers
```

---

## API reference

### `run_full_pipeline`

```python
from gwas_mr import run_full_pipeline

results = run_full_pipeline(
    disease_name: str,          # e.g. "Breast Cancer"
    biosample_type: str,        # e.g. "Breast Mammary Tissue", "PBMC", "Whole Blood"
    output_dir: str,            # MR results directory
    gwas_data_dir: str,         # GWAS downloads directory
    *,
    # All below are optional keyword arguments
    tsv_path: str = None,       # Custom GWAS Catalog TSV path
    eqtl_root: str = None,      # Custom eQTL directory
    gtex_biosample_csv: str = None, # Custom GTEx sample count CSV
    r_script_path: str = None,  # Custom mr_pipeline.R path
    plink_bin: str = None,      # Custom PLINK binary path
    plink_ref: str = None,      # Custom PLINK ref panel prefix
    use_llm: bool = True,       # Use OpenAI for disease normalization
    top_k: int = 1,             # Number of top GWAS studies per trait
    run_mr_analysis: bool = True,  # Set False to skip MR step
)
```

### Individual stages

```python
from gwas_mr import retrieve_gwas_data  # GWAS retrieval only
from gwas_mr import run_mr              # MR execution only
```

---

## Advanced configuration

### Custom paths (all optional)

```python
results = run_full_pipeline(
    disease_name="Type 2 Diabetes",
    biosample_type="Pancreas",
    output_dir="./results",
    gwas_data_dir="./gwas_cache",
    tsv_path="/path/to/GWAS-DATABASE-FTP.tsv",
    eqtl_root="/path/to/GTEx_eQTL_TISSUE_EXPRESSION",
    gtex_biosample_csv="/path/to/EQTL-SAMPLE-COUNT-GTEx Portal.csv",
    plink_bin=r"C:\plink\plink.exe",
    plink_ref=r"C:\ref_panels\g1000_eur",
)
```

### PLINK detection order

| Priority | Binary | Reference panel |
|---|---|---|
| 1 | `plink_bin` argument | `plink_ref` argument |
| 2 | `PLINK_BIN` env var | `PLINK_REF` env var |
| 3 | `plink` in system PATH | `gwas_mr_reference/plink/1000G_EUR_hg38` |
| 4 | `gwas_mr_reference/plink/plink.exe` or `plink` | — |

### Rscript detection order

| Priority | Location |
|---|---|
| 1 | `Rscript` in system PATH |
| 2 | `%R_HOME%\bin\Rscript.exe` (Windows) |
| 3 | `%R_HOME%\bin\x64\Rscript.exe` (Windows 64-bit) |

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `OPENAI_API_KEY not found` | Create a `.env` file with your key, or pass `use_llm=False` |
| `Rscript not found` | Add R to PATH, or set `R_HOME` env var to your R install directory |
| `PLINK status 5` | Ensure file paths don't have unquoted spaces; update to latest version |
| `Column 'se' not found` (R warning) | Already handled — the script auto-detects `standard_error` vs `se` |
| `Skipping clumping (<2 SNPs)` | Normal — gene has only 1 overlapping SNP, nothing to clump |
