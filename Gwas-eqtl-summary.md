Technical Briefing: GWAS and GTEx eQTL Mendelian Randomization Pipeline Navigator

Executive Summary

The GWAS and GTEx eQTL Mendelian Randomization (MR) Pipeline Navigator is a specialized script designed to automate the synthesis of Genome-Wide Association Study (GWAS) data and Genotype-Tissue Expression (GTEx) expression Quantitative Trait Loci (eQTL) data. Its primary function is to identify exact-matching European GWAS datasets for a specific disease or phenotype, map them to relevant GTEx tissue types, and generate ready-to-execute Mendelian Randomization commands.

The pipeline distinguishes itself through a multi-stage safety and normalization process, utilizing optional Large Language Model (LLM) integration for disease synonym identification and strict "Option A+" exact-match filtering to prevent the use of incorrect phenotype subtypes. By managing data acquisition through an "offline-first" cache system and automated downloading of harmonized summary statistics, the script significantly reduces the manual overhead required to prepare genetic epidemiology studies.


--------------------------------------------------------------------------------


Technical Workflow and Setup

The pipeline operates through a structured "pseudoprogram" sequence, beginning with environment configuration and ending with the generation of master reports.

Environmental Configuration

Before processing, the script establishes several foundational parameters:

* Environment Variables: Loads configurations (e.g., OPENAI_API_KEY) from a local .env file.
* Data Pathways: Defines paths for the GWAS catalog (GWAS-DATABASE-FTP.tsv), output directories, GTEx eQTL tissue files, and GTEx tissue sample-count CSVs.
* Operational Parameters: Sets limits for HTTP download retries, chunk sizes, and the number of top GWAS rows to consider (defined as TOP_K).

Workflow Overview

Step	Phase	Key Action
1	Input	User provides disease/phenotype and tissue names.
2	Normalization	(Optional) LLM-based synonym generation and name standardization.
3	Filtering	Identification of European-ancestry GWAS rows with available summary statistics.
4	Verification	"Option A+" safety filter ensures exact phenotype matching.
5	Tissue Mapping	Matching user input to GTEx tissue files using fuzzy or substring logic.
6	Retrieval	Checking local cache or downloading harmonized GWAS files.
7	Synthesis	Generating the run_pipeline.py command and comprehensive reports.


--------------------------------------------------------------------------------


Phenotype and Ancestry Optimization

The pipeline employs specific logic to ensure that only high-quality, relevant GWAS datasets are selected for the MR analysis.

Disease Name Normalization

If enabled via an API key, the script uses a GPT model to normalize the user's input. This produces a normalized_name and an array of synonyms[]. If the LLM is unavailable or the process fails, the script defaults to the raw user input for searching.

GWAS Table Processing

The script loads the GWAS FTP table and computes derived helper columns to refine selection:

* Ancestry Identification: The discoverySampleAncestry column is cleaned to retain only the primary ancestry.
* Sample Size (N): The script extracts the best-guess sample size from descriptive text, prioritizing "N=" strings or the maximum number found in the entry.
* European Cohort Filtering: A Boolean column is_EUR is created, flagging entries as TRUE only if the ancestry description confirms European descent.


--------------------------------------------------------------------------------


Data Integrity and Safety Filters

A critical component of the pipeline is the "Option A+" safety filter, which prevents the accidental inclusion of unrelated phenotypes or subtypes.

1. Broad Search: Initially, the script searches efoTraits for the normalized disease name or its synonyms.
2. Strict Enforcement: The script then enforces an exact match requirement (case-insensitive) between efoTraits and the normalized name.
3. Synonym Fallback: If no exact match is found for the primary name, the script pauses for five seconds before attempting an exact match against each synonym.
4. Failure Protocol: If no exact match is found after these steps, the script prints nearby available traits and terminates the process to maintain data precision.


--------------------------------------------------------------------------------


Tissue Mapping and eQTL Selection

The script maps user-requested tissues to specific GTEx files using a hierarchical matching logic:

* Matching Hierarchy: Exact match → Substring match → Fuzzy match.
* Default Behavior: If no match is found, the system defaults to "Whole Blood."
* Sample Size Retrieval: The script identifies both RNA-seq and Genotype sample sizes (N) from the GTEx tissue sample-count CSV to ensure accurate weighting in the MR command.


--------------------------------------------------------------------------------


Automated Data Acquisition and Caching

The pipeline prioritizes efficiency through an "offline-first" caching strategy for GWAS files.

* Cache Check: The script searches the output root for a sanitized folder matching the trait name. It detects existing harmonized files (.h.tsv.gz or .h.tsv).
* Automated Download: If no cache exists, the script sorts candidates by sample size N and attempts to download harmonized files (including meta.yaml and README) from the summary statistics URL.
* Resilience Features: The download process supports the Range header for resuming interrupted transfers and includes retry logic for broken connections.
* Format Conversion: To facilitate ease of use, downloaded .h.tsv.gz files are automatically unzipped to .h.tsv.


--------------------------------------------------------------------------------


Final Deliverables: MR Commands and Reports

The terminal stage of the pipeline synthesizes all gathered data into actionable outputs.

MR Command Generation

The script generates a ready-to-use Python command for the MR pipeline (run_pipeline.py). This command includes:

* Paths to the specific eQTL and GWAS files.
* Calculated sample sizes (n_eqtl and n_gwas), prioritizing Genotype N for eQTLs.
* Designated output directories for the results.

Comprehensive Reporting

The script produces two levels of documentation:

1. Dataset Excel Reports: Individual reports stored within each trait's folder, detailing accession numbers, ancestry, sample sizes, paths, and the specific run command used.
2. Master Reports: A global summary provided in both .tsv and .xlsx formats (MASTER_REPORT), aggregating all processed traits and their associated metadata.
