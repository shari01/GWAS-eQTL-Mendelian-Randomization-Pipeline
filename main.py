import argparse
from pathlib import Path

import pandas as pd

from gwas_mr import run_full_pipeline


REPO_ROOT = Path(__file__).resolve().parent
ASSOC_REQUIRED_COLUMNS = {"DISEASE/TRAIT", "STUDY", "REPORTED GENE(S)"}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the GWAS-eQTL MR pipeline from this repository."
    )
    parser.add_argument(
        "deg_file",
        nargs="?",
        help="Path to a DEG/gene CSV/TSV/XLSX file for the GWAS association step.",
    )
    parser.add_argument(
        "--disease",
        default="Acute Myeloid Leukemia",
        help="Disease name used for GWAS retrieval and MR.",
    )
    parser.add_argument(
        "--biosample",
        default="Bone Marrow",
        help="GTEx biosample/tissue to pair with the disease.",
    )
    parser.add_argument(
        "--output-dir",
        default="pipeline_output_test",
        help="Directory for MR and association outputs.",
    )
    parser.add_argument(
        "--gwas-data-dir",
        default="GWAS_Data_Output",
        help="Directory for cached/downloaded GWAS datasets.",
    )
    parser.add_argument(
        "--disable-llm",
        action="store_true",
        help="Skip OpenAI-based disease normalization.",
    )
    parser.add_argument(
        "--skip-preflight",
        action="store_true",
        help="Skip dependency checks before running.",
    )
    return parser


def has_valid_assoc_index(path: Path) -> bool:
    """Return True when the association index looks compatible."""
    if not path.is_file():
        return False
    try:
        header = pd.read_csv(path, sep="\t", nrows=0, low_memory=False)
    except Exception:
        return False
    return ASSOC_REQUIRED_COLUMNS.issubset(set(header.columns))


def main() -> None:
    args = build_parser().parse_args()

    deg_path = Path(args.deg_file).expanduser().resolve() if args.deg_file else None
    if deg_path and not deg_path.is_file():
        raise FileNotFoundError(f"DEG file not found: {deg_path}")

    assoc_index = REPO_ROOT / "gwas-latest-all-assoc.tsv"
    run_association = bool(deg_path and has_valid_assoc_index(assoc_index))

    if deg_path and not assoc_index.is_file():
        print(
            "[WARN] gwas-latest-all-assoc.tsv not found in the repo root. "
            "Running MR only and skipping the GWAS association step.",
            flush=True,
        )
    elif deg_path and assoc_index.is_file() and not run_association:
        print(
            "[WARN] gwas-latest-all-assoc.tsv exists but does not match the "
            "expected GWAS association format. Running MR only.",
            flush=True,
        )

    print(f"[RUN] disease={args.disease} | biosample={args.biosample}", flush=True)
    print(f"[RUN] deg_file={deg_path if deg_path else 'none'}", flush=True)

    run_full_pipeline(
        disease_name=args.disease,
        biosample_type=args.biosample,
        input_genes=deg_path,
        output_dir=str(REPO_ROOT / args.output_dir),
        gwas_data_dir=str(REPO_ROOT / args.gwas_data_dir),
        tsv_path=str(REPO_ROOT / "GWAS-DATABASE-FTP.tsv"),
        eqtl_root=str(REPO_ROOT / "GTEx_eQTL_TISSUE_EXPRESSION"),
        gtex_biosample_csv=str(
            REPO_ROOT / "GTEx_EQTL_SAMPLE_COUNTS" / "EQTL-SAMPLE-COUNT-GTEx Portal.csv"
        ),
        plink_bin=str(REPO_ROOT / "plink" / "plink.exe"),
        plink_ref=str(REPO_ROOT / "plink" / "1000G_EUR_hg38"),
        gwas_assoc_index=str(assoc_index),
        use_llm=not args.disable_llm,
        run_association=run_association,
        skip_preflight=args.skip_preflight,
    )


if __name__ == "__main__":
    main()


