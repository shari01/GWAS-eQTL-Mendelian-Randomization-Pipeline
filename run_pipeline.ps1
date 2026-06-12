<#
.SYNOPSIS
    Run the GWAS-eQTL Mendelian Randomization pipeline.

.DESCRIPTION
    Wraps main.py via uv.  All paths can be absolute or relative to this script.

    -Biosample accepts TWO forms:
      (A) AUTO — any tissue/fluid description; LLM maps it to the closest GTEx tissue.
      (B) MANUAL — an exact filename from GTEx_eQTL_TISSUE_EXPRESSION\; the file is
          used directly and the tissue name is parsed from the filename.

.EXAMPLE  AUTO AI tissue resolution (LLM maps "urine" → kidney cortex at temperature=0)
    .\run_pipeline.ps1 "D:\AyassBio_Workspace_Downloads\sabhat-dataset-tr\lupus-degs.csv" `
        -Disease "diabetes" `
        -Biosample "urine"

.EXAMPLE  AUTO — other fluid/cell type examples (all resolved by AI)
    .\run_pipeline.ps1 "degs.csv" -Disease "Alzheimer's disease" -Biosample "CSF"
    .\run_pipeline.ps1 "degs.csv" -Disease "lupus"               -Biosample "PBMC"
    .\run_pipeline.ps1 "degs.csv" -Disease "Crohn's disease"     -Biosample "stool"

.EXAMPLE  MANUAL — pass an exact GTEx eQTL filename (skips AI resolution entirely)
    .\run_pipeline.ps1 "D:\AyassBio_Workspace_Downloads\sabhat-dataset-tr\lupus-degs.csv" `
        -Disease "diabetes" `
        -Biosample "Kidney_Cortex.v10.eGenes.txt"

.EXAMPLE  MANUAL — other tissue file examples
    .\run_pipeline.ps1 "degs.csv" -Disease "Alzheimer's disease" -Biosample "Brain_Cortex.v10.eGenes.txt"
    .\run_pipeline.ps1 "degs.csv" -Disease "lupus"               -Biosample "Whole_Blood.v10.eGenes.txt"
    .\run_pipeline.ps1 "degs.csv" -Disease "type 2 diabetes"     -Biosample "Pancreas.v10.eGenes.txt"

.EXAMPLE  MR only (no DEG/gene file), disable LLM for speed
    .\run_pipeline.ps1 -Disease "coronary artery disease" -Biosample "Artery_Coronary.v10.eGenes.txt" -DisableLlm

.EXAMPLE  Skip preflight dependency checks
    .\run_pipeline.ps1 "degs.csv" -Disease "breast cancer" -Biosample "Breast_Mammary_Tissue.v10.eGenes.txt" -SkipPreflight

.NOTES
    Available GTEx eQTL tissue files (pass any of these as -Biosample for manual mode):
      Adipose_Subcutaneous.v10.eGenes.txt
      Adipose_Visceral_Omentum.v10.eGenes.txt
      Adrenal_Gland.v10.eGenes.txt
      Artery_Aorta.v10.eGenes.txt
      Artery_Coronary.v10.eGenes.txt
      Artery_Tibial.v10.eGenes.txt
      Bladder.v10.eGenes.txt
      Brain_Amygdala.v10.eGenes.txt
      Brain_Anterior_cingulate_cortex_BA24.v10.eGenes.txt
      Brain_Caudate_basal_ganglia.v10.eGenes.txt
      Brain_Cerebellar_Hemisphere.v10.eGenes.txt
      Brain_Cerebellum.v10.eGenes.txt
      Brain_Cortex.v10.eGenes.txt
      Brain_Frontal_Cortex_BA9.v10.eGenes.txt
      Brain_Hippocampus.v10.eGenes.txt
      Brain_Hypothalamus.v10.eGenes.txt
      Brain_Nucleus_accumbens_basal_ganglia.v10.eGenes.txt
      Brain_Putamen_basal_ganglia.v10.eGenes.txt
      Brain_Spinal_cord_cervical_c-1.v10.eGenes.txt
      Brain_Substantia_nigra.v10.eGenes.txt
      Breast_Mammary_Tissue.v10.eGenes.txt
      Cells_Cultured_fibroblasts.v10.eGenes.txt
      Cells_EBV-transformed_lymphocytes.v10.eGenes.txt
      Colon_Sigmoid.v10.eGenes.txt
      Colon_Transverse.v10.eGenes.txt
      Esophagus_Gastroesophageal_Junction.v10.eGenes.txt
      Esophagus_Mucosa.v10.eGenes.txt
      Esophagus_Muscularis.v10.eGenes.txt
      Heart_Atrial_Appendage.v10.eGenes.txt
      Heart_Left_Ventricle.v10.eGenes.txt
      Kidney_Cortex.v10.eGenes.txt
      Liver.v10.eGenes.txt
      Lung.v10.eGenes.txt
      Minor_Salivary_Gland.v10.eGenes.txt
      Muscle_Skeletal.v10.eGenes.txt
      Nerve_Tibial.v10.eGenes.txt
      Ovary.v10.eGenes.txt
      Pancreas.v10.eGenes.txt
      Pituitary.v10.eGenes.txt
      Prostate.v10.eGenes.txt
      Skin_Not_Sun_Exposed_Suprapubic.v10.eGenes.txt
      Skin_Sun_Exposed_Lower_leg.v10.eGenes.txt
      Small_Intestine_Terminal_Ileum.v10.eGenes.txt
      Spleen.v10.eGenes.txt
      Stomach.v10.eGenes.txt
      Testis.v10.eGenes.txt
      Thyroid.v10.eGenes.txt
      Uterus.v10.eGenes.txt
      Vagina.v10.eGenes.txt
      Whole_Blood.v10.eGenes.txt
#>
param(
    [Parameter(Position = 0)]
    [string]$DegFile,

    [string]$Disease = "Acute Myeloid Leukemia",

    # Either a tissue/fluid description (AI resolves it) OR an exact GTEx eQTL filename.
    # AUTO:   -Biosample "urine"                              (LLM maps to kidney cortex)
    # MANUAL: -Biosample "Kidney_Cortex.v10.eGenes.txt"      (uses that file directly)
    [string]$Biosample = "Bone Marrow",

    [string]$OutputDir = "pipeline_output_test",

    [string]$GwasDataDir = "GWAS_Data_Output",

    [switch]$DisableLlm,

    [switch]$SkipPreflight
)

# Clear any shell-activated venv so uv uses its own managed environment without warning.
$env:VIRTUAL_ENV = $null

$argsList = @(
    "run",
    "--python", "3.12",
    "--with-requirements", "requirements.txt",
    "--env-file", ".env",
    "python",
    "main.py"
)

if ($DegFile) {
    $argsList += $DegFile
}

$argsList += @(
    "--disease", $Disease,
    "--biosample", $Biosample,
    "--output-dir", $OutputDir,
    "--gwas-data-dir", $GwasDataDir
)

if ($DisableLlm) {
    $argsList += "--disable-llm"
}

if ($SkipPreflight) {
    $argsList += "--skip-preflight"
}

& uv @argsList
exit $LASTEXITCODE
