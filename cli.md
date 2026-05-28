# CLI Run Commands

## Main Run Command

### PowerShell

```powershell
.\run_pipeline.ps1 "C:\Users\shahr\Downloads\GWAS-eQTL Mendelian Randomization Pipeline\Chronic_Kidney_Disease_100_genes (1).csv" -Disease "Chronic Kidney Disease" -Biosample "Kidney Tissue"
```

### Command Prompt

```cmd
run_pipeline.cmd "C:\Users\shahr\Downloads\GWAS-eQTL Mendelian Randomization Pipeline\Chronic_Kidney_Disease_100_genes (1).csv" --disease "Chronic Kidney Disease" --biosample "Kidney Tissue"
```

## Help Commands

### PowerShell

```powershell
.\run_pipeline.ps1 --help
```

### Command Prompt

```cmd
run_pipeline.cmd --help
```

## Direct Python Command

```powershell
uv run --python 3.12 --with-requirements requirements.txt --env-file .env python main.py "C:\Users\shahr\Downloads\GWAS-eQTL Mendelian Randomization Pipeline\Chronic_Kidney_Disease_100_genes (1).csv" --disease "Chronic Kidney Disease" --biosample "Kidney Tissue"
```

## What Each Input Means

- First file path:
  Your DEG file path.
- `-Disease` or `--disease`:
  Disease name for GWAS retrieval and MR.
- `-Biosample` or `--biosample`:
  GTEx tissue/biosample name.

## Output Locations

- `GWAS_Data_Output/`
  GWAS retrieval and cache output
- `pipeline_output_test/`
  MR results and reports

## Notes

- `main.py` is the real entry point.
- `run_pipeline.ps1` is the easiest command for PowerShell.
- `run_pipeline.cmd` is the easiest command for Command Prompt.
- If `gwas-latest-all-assoc.tsv` is present and valid, the GWAS association step also runs.
