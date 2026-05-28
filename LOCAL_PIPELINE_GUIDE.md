# Local Pipeline Guide

## What This Repo Runs

This repository runs an end-to-end GWAS + GTEx eQTL + Mendelian Randomization workflow.

When you launch [main.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/main.py:1), the pipeline does:

1. Reads your disease name and biosample/tissue.
2. Uses `GWAS-DATABASE-FTP.tsv` to find the best matching GWAS study.
3. Uses `GTEx_eQTL_TISSUE_EXPRESSION/` and `GTEx_EQTL_SAMPLE_COUNTS/` for tissue-specific eQTL data.
4. Runs the R MR pipeline in [gwas_mr/mr_pipeline.R](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/mr_pipeline.R:1).
5. Uses `plink/1000G_EUR_hg38*` for LD-aware steps.
6. If `gwas-latest-all-assoc.tsv` is present and valid, also runs the GWAS association sub-pipeline using your DEG/gene file.

## Main Entry Points

- [main.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/main.py:1)
  Main Python entry point. Accepts DEG file path plus disease and biosample arguments.

- [run_pipeline.ps1](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/run_pipeline.ps1:1)
  PowerShell wrapper for easy local runs.

- [run_pipeline.cmd](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/run_pipeline.cmd:1)
  Windows Command Prompt wrapper for easy local runs.

## Key Code Files

- [gwas_mr/defaults.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/defaults.py:1)
  Resolves reference-data paths from this repo root.

- [gwas_mr/pipeline.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/pipeline.py:290)
  Main orchestrator. Runs preflight, GWAS retrieval, MR, and optional association.

- [gwas_mr/retrieval.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/retrieval.py:652)
  GWAS dataset lookup, tissue matching, GWAS download/cache logic.

- [gwas_mr/mr_runner.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/mr_runner.py:129)
  Calls `Rscript` and passes PLINK paths to the R backend.

- [gwas_mr/report.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/report.py:394)
  Builds the HTML MR summary report.

- [gwas_mr/gwas_association/pipeline.py](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/gwas_mr/gwas_association/pipeline.py:73)
  Optional association step using your DEG/gene file plus `gwas-latest-all-assoc.tsv`.

## Required Local Files

The current code expects these large files and folders in the repository root:

- `GWAS-DATABASE-FTP.tsv`
- `gwas-latest-all-assoc.tsv`
- `GTEx_eQTL_TISSUE_EXPRESSION/`
- `GTEx_EQTL_SAMPLE_COUNTS/`
- `plink/`
- `.env`

Inside `plink/`, the working LD reference set is:

- `1000G_EUR_hg38.bed`
- `1000G_EUR_hg38.bim`
- `1000G_EUR_hg38.fam`

## How To Run

### PowerShell

```powershell
.\run_pipeline.ps1 "C:\path\to\your_deg_file.csv" -Disease "Chronic Kidney Disease" -Biosample "Kidney Tissue"
```

### Command Prompt

```cmd
run_pipeline.cmd "C:\path\to\your_deg_file.csv" --disease "Chronic Kidney Disease" --biosample "Kidney Tissue"
```

### Direct Python command

```powershell
uv run --python 3.12 --with-requirements requirements.txt --env-file .env python main.py "C:\path\to\your_deg_file.csv" --disease "Chronic Kidney Disease" --biosample "Kidney Tissue"
```

## Example For This Repo

```powershell
.\run_pipeline.ps1 "C:\Users\shahr\Downloads\GWAS-eQTL Mendelian Randomization Pipeline\Chronic_Kidney_Disease_100_genes (1).csv" -Disease "Chronic Kidney Disease" -Biosample "Kidney Tissue"
```

## What Inputs Mean

- DEG file
  A CSV, TSV, TXT, XLSX, or XLS file with a gene column such as `gene`, `gene_symbol`, or `symbol`.

- Disease
  Used for GWAS study lookup and report labels.

- Biosample
  Used to select the GTEx tissue eQTL file. The code will fuzzy-match this to GTEx tissue names.

## Outputs

- `GWAS_Data_Output/`
  Cached or downloaded GWAS study folders plus a master report.

- `pipeline_output_test/<trait>/`
  MR results, diagnostics, summary tables, and `MR_PIPELINE_REPORT.html`.

- `pipeline_output_test/<disease>/gwas_association/`
  Association CSVs and HTML graphs, when the association index is valid and present.

## Current Validation Status

The repo has already been validated for:

- Python dependency loading through `uv`
- `.env` loading
- `Rscript` availability
- required R packages
- PLINK executable
- PLINK LD reference prefix
- `GWAS-DATABASE-FTP.tsv`
- `gwas-latest-all-assoc.tsv`

## Quick Troubleshooting

- If the DEG file path is wrong, `main.py` stops immediately with `DEG file not found`.
- If `gwas-latest-all-assoc.tsv` is missing or invalid, MR still runs but the association branch is skipped.
- If `OPENAI_API_KEY` is missing from `.env`, disable LLM mode or add the key back.
- Large outputs are written locally and are ignored by git via [.gitignore](/C:/Users/shahr/Downloads/GWAS-eQTL%20Mendelian%20Randomization%20Pipeline/.gitignore:1).
