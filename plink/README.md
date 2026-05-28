# PLINK Setup for MR Pipeline

The MR pipeline uses **PLINK 1.9** for LD clumping and LD matrix
computation. Without it, the pipeline still runs but skips those steps.

## 1. Install PLINK 1.9

Download from: https://www.cog-genomics.org/plink/

- **Windows**: download `plink_win64.zip`, extract `plink.exe` into this
  directory (`gwas_ret/data/plink/plink.exe`)
- **Linux**: download the Linux binary, extract `plink` into this directory
- **macOS**: download the macOS binary, extract `plink` into this directory

Or install system-wide and add to your PATH — the pipeline auto-detects it.

## 2. Download the LD Reference Panel (1000 Genomes EUR, hg38)

The pipeline needs a 1000 Genomes European reference panel in PLINK
binary format (`.bed`, `.bim`, `.fam`) built on **GRCh38/hg38**.

Place the three files in this directory with the prefix `1000G_EUR_hg38`:

```
gwas_ret/data/plink/
├── plink.exe            (or plink on Linux/macOS)
├── 1000G_EUR_hg38.bed
├── 1000G_EUR_hg38.bim
└── 1000G_EUR_hg38.fam
```

### Where to get the reference panel

- **PLINK 2 resources page**: https://www.cog-genomics.org/plink/2.0/resources
- **1000 Genomes on GRCh38**: filter for EUR samples and convert to PLINK
  format, or use pre-built files from the link above.

## 3. Alternative: Environment Variables

Instead of placing files here, you can set environment variables:

```
PLINK_BIN=C:\path\to\plink.exe
PLINK_REF=C:\path\to\1000G_EUR_hg38
```

Or pass them directly in Python:

```python
from gwas_ret import run_full_pipeline

results = run_full_pipeline(
    disease_name="Breast Cancer",
    tissue_name="Whole Blood",
    output_dir="./output",
    plink_bin=r"C:\plink\plink.exe",
    plink_ref=r"D:\ref_panels\1000G_EUR_hg38",
)
```

## Detection order

The pipeline searches for PLINK in this order:

| Priority | PLINK binary | Reference panel |
|---|---|---|
| 1 | `plink_bin` function argument | `plink_ref` function argument |
| 2 | `PLINK_BIN` env var | `PLINK_REF` env var |
| 3 | `plink` in system PATH | `gwas_ret/data/plink/1000G_EUR_hg38` |
| 4 | `gwas_ret/data/plink/plink.exe` | — |
| 5 | *(skip — MR runs without LD)* | — |
