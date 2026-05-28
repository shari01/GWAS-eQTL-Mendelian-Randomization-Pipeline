param(
    [Parameter(Position = 0)]
    [string]$DegFile,

    [string]$Disease = "Acute Myeloid Leukemia",

    [string]$Biosample = "Bone Marrow",

    [string]$OutputDir = "pipeline_output_test",

    [string]$GwasDataDir = "GWAS_Data_Output",

    [switch]$DisableLlm,

    [switch]$SkipPreflight
)

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
