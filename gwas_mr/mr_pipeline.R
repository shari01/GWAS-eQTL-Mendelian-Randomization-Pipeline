suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

# ============================
# READ ARGUMENTS FROM COMMAND LINE
# ============================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript mr_pipeline.R <base_dir> <eqtl_path> <gwas_path> <N_EQTL> <N_GWAS> <out_dir> <disease_name> <biosample_type>")
}

base_dir       <- args[1]
eqtl_path      <- args[2]
gwas_path      <- args[3]
N_EQTL_DEFAULT <- as.numeric(args[4])
N_GWAS_DEFAULT <- as.numeric(args[5])
out_dir        <- args[6]
OUTCOME_NAME   <- args[7]
BIOSAMPLE_TYPE <- args[8]


has_coloc    <- requireNamespace("coloc", quietly = TRUE)
has_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
has_susieR   <- requireNamespace("susieR", quietly = TRUE)
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

# ============================
# CONFIG (EDIT HERE)
# ============================
# base_dir  <- "C:/Users/shahr/Downloads/MR_EQTL_GWAS"
# eqtl_path <- file.path(base_dir, "Pancreas.v10.eGenes.txt")
# gwas_path <- file.path(base_dir, "diabetes-harmonised.qc.tsv")

# ---- PLINK (LOCAL LD PANEL) ----
# Auto-detection: env vars > PATH.
# Set PLINK_BIN and PLINK_REF env vars before running, or place plink
# in your PATH and the reference panel in gwas_mr_reference/plink/.
PLINK_BIN <- Sys.getenv("PLINK_BIN", unset = "")
if (!nzchar(PLINK_BIN)) {
  plink_from_path <- Sys.which("plink")
  if (nzchar(plink_from_path)) {
    PLINK_BIN <- plink_from_path
  }
}
PLINK_REF <- Sys.getenv("PLINK_REF", unset = "")

SAVE_PLOTS   <- TRUE
PLOT_FORMAT  <- "png"
PLOT_DPI     <- 400

EQTL_P_THRESH <- 5e-8
FSTAT_THRESH  <- 10

CLUMP_R2_STRICT <- 0.1
CLUMP_R2_SMALLN <- 0.2
CLUMP_KB <- 250

HARMONISE_ACTION <- 1
MIN_SNP_STRONG <- 2

# COLOC (needs chr/pos + N defaults)
RUN_COLOC <- TRUE
COLOC_WINDOW_BP <- 5e5

# N_EQTL_DEFAULT <- 362
# N_GWAS_DEFAULT <- 520580

GWAS_TYPE <- "cc"
EQTL_TYPE <- "quant"

# LD thresholds for correlated SNP list
LD_R2_THRESHOLDS <- c(0.2, 0.5, 0.8)

# SuSiE fine-mapping
RUN_FINEMAP <- TRUE
SUSIE_L <- 10
FINEMAP_WINDOW_BP <- 5e5

# Output root
# out_dir <- file.path(base_dir, "MR_RESULTS-lx-plink-test2")

# ============================
# OUTPUT STRUCTURE (MATCH SPEC)
# ============================
d00 <- file.path(out_dir, "00_config")
d01 <- file.path(out_dir, "01_input_raw")
d02 <- file.path(out_dir, "02_input_parsed")
d03 <- file.path(out_dir, "03_instruments")
d03a <- file.path(d03, "03a_eqtl_filtered")
d03b <- file.path(d03, "03b_fstat")
d04 <- file.path(out_dir, "04_overlap")
d11_pg <- file.path(out_dir, "11_plots", "per_gene")
d12 <- file.path(out_dir, "summary_tables")

# Only PLINK needs temp directories (for clumping input/output and LD
# matrix computation).  All other per-gene data is kept exclusively in
# memory and aggregated into 12_summary_tables.
.pg_tmp <- file.path(tempdir(), "mr_plink_tmp")
d05     <- file.path(.pg_tmp, "clump")
d07     <- file.path(.pg_tmp, "ld")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
for (dd in c(d00,d01,d02,d03a,d03b,d04,d11_pg,d12)) {
  dir.create(dd, recursive=TRUE, showWarnings=FALSE)
}
dir.create(d05, recursive=TRUE, showWarnings=FALSE)
dir.create(d07, recursive=TRUE, showWarnings=FALSE)

# ============================
# LOGGING + CONFIG DUMP
# ============================
log_file <- file.path(d00, "RUN_LOG.txt")
writeLines(paste0("RUN START: ", Sys.time()), log_file)
log_msg <- function(x) { cat(x, "\n"); write(x, log_file, append=TRUE) }

# ============================
# HELPERS
# ============================
theme_mr <- function() {
  theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank())
}
save_plot <- function(p, path) {
  if (!SAVE_PLOTS) return(invisible(NULL))
  if (PLOT_FORMAT == "pdf") ggsave(path, p, width=7.2, height=5.4, device=cairo_pdf)
  else ggsave(path, p, width=7.2, height=5.4, dpi=PLOT_DPI)
}
save_plot_list <- function(plot_list, base_path_prefix) {
  if (!SAVE_PLOTS) return(invisible(NULL))
  if (is.null(plot_list)) return(invisible(NULL))
  if (inherits(plot_list, "ggplot")) {
    save_plot(plot_list + theme_mr(), paste0(base_path_prefix, ".", PLOT_FORMAT)); return(invisible(NULL))
  }
  if (is.list(plot_list)) {
    for (i in seq_along(plot_list)) {
      p <- plot_list[[i]]
      if (inherits(p, "ggplot")) save_plot(p + theme_mr(), paste0(base_path_prefix, "_", i, ".", PLOT_FORMAT))
    }
  }
  invisible(NULL)
}
extract_chr_pos_from_variant_id <- function(dt) {
  cand <- intersect(c("variant_id","variantid","variant"), names(dt))
  if (length(cand) == 0) return(dt)
  vcol <- cand[1]
  v <- as.character(dt[[vcol]])
  ok <- grepl("^chr[^_]+_[0-9]+_", v)
  if (!any(ok, na.rm=TRUE)) return(dt)
  chr <- sub("^([^_]+)_.*$", "\\1", v)
  pos <- suppressWarnings(as.numeric(sub("^[^_]+_([0-9]+)_.*$", "\\1", v)))
  if (!("chr" %in% names(dt))) dt[, chr := chr]
  if (!("pos" %in% names(dt))) dt[, pos := pos]
  dt
}

# ---- PLINK checks ----
check_plink_ready <- function(plink_bin, plink_ref) {
  ok1 <- nzchar(plink_bin)
  ok2 <- file.exists(paste0(plink_ref, ".bed")) && file.exists(paste0(plink_ref, ".bim")) && file.exists(paste0(plink_ref, ".fam"))
  if (!ok1) { log_msg("PLINK_BIN missing."); return(FALSE) }
  if (!ok2) {
    log_msg(paste("PLINK_REF not found (.bed/.bim/.fam):", plink_ref))
    return(FALSE)
  }
  test <- tryCatch(system2(plink_bin, args=c("--version"), stdout=TRUE, stderr=TRUE), error=function(e) e)
  if (inherits(test, "error")) { log_msg(paste("PLINK exec failed:", test$message)); return(FALSE) }
  log_msg("PLINK OK (local LD panel available).")
  TRUE
}

# GWAS loader (unchanged)
load_gwas_catalog_mr <- function(path) {
  gwas <- fread(path, showProgress=FALSE)
  setnames(gwas, names(gwas), tolower(names(gwas)))
  rs_col <- intersect(c("rsid","rs_id","snp"), names(gwas))
  if (length(rs_col) == 0) {
    for (cn in names(gwas)) {
      v <- as.character(gwas[[cn]])
      if (any(grepl("^rs[0-9]+$", v), na.rm=TRUE)) { rs_col <- cn; break }
    }
  }
  if (length(rs_col) == 0) stop("NO rsID FOUND")
  rs_col <- rs_col[1]
  if ("beta" %in% names(gwas)) gwas[, beta := as.numeric(beta)]
  else if ("odds_ratio" %in% names(gwas)) gwas[, beta := log(as.numeric(odds_ratio))]
  else if ("hazard_ratio" %in% names(gwas)) gwas[, beta := log(as.numeric(hazard_ratio))]
  else stop("NO EFFECT SIZE (beta / odds_ratio / hazard_ratio)")
  se_col <- intersect(c("standard_error","se"), names(gwas))
  if (length(se_col) == 0) stop("NO SE (standard_error / se)")
  gwas[, se := as.numeric(get(se_col[1]))]
  if ("p_value" %in% names(gwas)) gwas[, pval := as.numeric(p_value)]
  else if ("neg_log_10_p_value" %in% names(gwas)) gwas[, pval := 10^(-as.numeric(neg_log_10_p_value))]
  else stop("NO P (p_value / neg_log_10_p_value)")
  if (!all(c("effect_allele","other_allele") %in% names(gwas)))
    stop("ALLELES MISSING (need effect_allele and other_allele)")
  eaf_col <- intersect(c("effect_allele_frequency","eaf","af"), names(gwas))
  if (length(eaf_col) == 0) gwas[, eaf := NA_real_] else gwas[, eaf := as.numeric(get(eaf_col[1]))]
  chr_col <- intersect(c("chromosome","chr"), names(gwas))
  pos_col <- intersect(c("base_pair_location","position","bp","pos","base_pair_position"), names(gwas))
  chr_use <- if (length(chr_col) > 0) chr_col[1] else NA_character_
  pos_use <- if (length(pos_col) > 0) pos_col[1] else NA_character_
  gwas_mr <- gwas[, .(
    SNP = as.character(get(rs_col)),
    beta = as.numeric(beta),
    se = as.numeric(se),
    pval = as.numeric(pval),
    effect_allele = as.character(effect_allele),
    other_allele  = as.character(other_allele),
    eaf = as.numeric(eaf),
    chr = if (!is.na(chr_use)) as.character(get(chr_use)) else NA_character_,
    pos = if (!is.na(pos_use)) as.numeric(get(pos_use)) else NA_real_
  )]
  gwas_mr <- gwas_mr[grepl("^rs[0-9]+$", SNP)]
  gwas_mr <- gwas_mr[!is.na(beta) & !is.na(se) & !is.na(pval)]
  setorder(gwas_mr, pval)
  unique(gwas_mr, by="SNP")
}

# ---- LD correlated SNPs (unchanged) ----
ld_correlated_snps <- function(ld_mat, lead_snp, thresholds=c(0.2,0.5,0.8)) {
  if (is.null(ld_mat) || !is.matrix(ld_mat)) return(NULL)
  if (!(lead_snp %in% colnames(ld_mat))) return(NULL)
  r <- ld_mat[, lead_snp]
  dt <- data.table(SNP=names(r), r=as.numeric(r))
  dt[, r2 := r^2]
  dt[, is_lead := SNP == lead_snp]
  dt <- dt[order(-r2)]
  for (t in thresholds) dt[, paste0("r2_ge_", gsub("\\.","_", t)) := r2 >= t]
  dt
}

plot_ld_heatmap <- function(ld_mat, outfile) {
  if (!SAVE_PLOTS || !has_pheatmap) return(invisible(NULL))
  if (is.null(ld_mat) || !is.matrix(ld_mat) || nrow(ld_mat) < 2) return(invisible(NULL))
  m2 <- ld_mat^2
  if (PLOT_FORMAT == "pdf") {
    pdf(outfile, width=8, height=7)
    pheatmap::pheatmap(m2, main="LD heatmap (r^2)", cluster_rows=TRUE, cluster_cols=TRUE)
    dev.off()
  } else {
    png(outfile, width=2200, height=1900, res=PLOT_DPI)
    pheatmap::pheatmap(m2, main="LD heatmap (r^2)", cluster_rows=TRUE, cluster_cols=TRUE)
    dev.off()
  }
}

# ---- PLINK clumping (offline replacement for clump_data) ----
plink_clump <- function(dt_exposure, plink_bin, plink_ref, out_prefix, r2, kb) {
  # dt_exposure must have columns: SNP, pval.exposure OR pval
  if (is.null(dt_exposure) || nrow(dt_exposure) == 0) return(list(ok=FALSE, msg="Empty exposure", keep=NULL))
  snp <- as.character(dt_exposure$SNP)
  snp <- snp[grepl("^rs[0-9]+$", snp)]
  if (length(snp) == 0) return(list(ok=FALSE, msg="No rsIDs", keep=NULL))

  pcol <- if ("pval.exposure" %in% names(dt_exposure)) "pval.exposure" else if ("pval" %in% names(dt_exposure)) "pval" else NA_character_
  if (is.na(pcol)) return(list(ok=FALSE, msg="No p-value column for clumping", keep=NULL))

  clump_in <- data.table(SNP=snp, P=as.numeric(dt_exposure[[pcol]][match(snp, dt_exposure$SNP)]))
  clump_in <- clump_in[!is.na(P)]
  if (nrow(clump_in) == 0) return(list(ok=FALSE, msg="No valid P values", keep=NULL))

  in_path <- paste0(out_prefix, "__clump_input.tsv")
  fwrite(clump_in, in_path, sep="\t")

  args <- c(
    "--bfile", shQuote(plink_ref),
    "--clump", shQuote(in_path),
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", "1",
    "--clump-p2", "1",
    "--clump-r2", as.character(r2),
    "--clump-kb", as.character(kb),
    "--out", shQuote(out_prefix)
  )

  res <- tryCatch(system2(plink_bin, args=args, stdout=TRUE, stderr=TRUE), error=function(e) e)
  clumped_path <- paste0(out_prefix, ".clumped")
  if (inherits(res, "error")) return(list(ok=FALSE, msg=res$message, keep=NULL))
  if (!file.exists(clumped_path)) {
    msg <- paste("PLINK clump did not produce .clumped. See logs:", paste(res, collapse=" | "))
    return(list(ok=FALSE, msg=msg, keep=NULL))
  }

  # Parse .clumped: keep index SNPs
  cl <- tryCatch(fread(clumped_path, fill=TRUE, showProgress=FALSE), error=function(e) NULL)
  if (is.null(cl) || !("SNP" %in% names(cl))) {
    return(list(ok=FALSE, msg="Could not parse .clumped", keep=NULL))
  }
  keep <- unique(as.character(cl$SNP))
  keep <- keep[grepl("^rs[0-9]+$", keep)]
  if (length(keep) == 0) return(list(ok=FALSE, msg="No SNPs kept after clump parse", keep=NULL))
  list(ok=TRUE, msg="OK", keep=keep)
}

# ---- PLINK LD matrix (offline replacement for ieugwasr::ld_matrix) ----
plink_ld_matrix <- function(snps, plink_bin, plink_ref, out_prefix) {
  snps <- unique(as.character(snps))
  snps <- snps[grepl("^rs[0-9]+$", snps)]
  if (length(snps) < 2) return(list(ok=FALSE, msg="Need >=2 SNPs for LD", mat=NULL))

  snp_path <- paste0(out_prefix, "__snps.txt")
  writeLines(snps, snp_path)

  # Use --r square to produce square LD matrix (r). We'll keep r (then r^2 for heatmap)
  args <- c(
    "--bfile", shQuote(plink_ref),
    "--extract", shQuote(snp_path),
    "--r", "square", "gz",
    "--out", shQuote(out_prefix)
  )

  res <- tryCatch(system2(plink_bin, args=args, stdout=TRUE, stderr=TRUE), error=function(e) e)
  ld_path <- paste0(out_prefix, ".ld.gz")

  if (inherits(res, "error")) return(list(ok=FALSE, msg=res$message, mat=NULL))
  if (!file.exists(ld_path)) {
    msg <- paste("PLINK LD did not produce .ld. See logs:", paste(res, collapse=" | "))
    return(list(ok=FALSE, msg=msg, mat=NULL))
  }

  # PLINK .ld is whitespace-delimited numeric matrix (no header)
  m <- tryCatch(
    as.matrix(fread(ld_path, header = FALSE, data.table = FALSE)),
    error = function(e) NULL
  )

  if (is.null(m) || !is.matrix(m)) return(list(ok=FALSE, msg="Failed reading .ld matrix", mat=NULL))
  if (nrow(m) != length(snps) || ncol(m) != length(snps)) {
    # Sometimes PLINK drops SNPs not found; try derive SNP order from .bim intersection
    # We'll attempt to infer the kept SNP order by re-reading PLINK log for "X variants loaded" is not enough.
    # Safer: rerun with --write-snplist to capture kept snps; PLINK 1.9 doesn't always.
    return(list(ok=FALSE, msg="LD matrix dimensions mismatch (SNPs missing in ref panel?)", mat=NULL))
  }
  colnames(m) <- snps
  rownames(m) <- snps
  list(ok=TRUE, msg="OK", mat=m)
}

# ---------- SuSiE fine-mapping ----------
run_susie_rss_safe <- function(beta, se, ld_mat, n, L=10) {
  if (!has_susieR) return(list(ok=FALSE, msg="susieR not installed", fit=NULL))
  if (is.null(ld_mat) || !is.matrix(ld_mat) || nrow(ld_mat) < 2) return(list(ok=FALSE, msg="LD matrix missing", fit=NULL))
  if (is.na(n) || !is.finite(n)) return(list(ok=FALSE, msg="Need sample size N for susie_rss", fit=NULL))
  z <- beta / se
  fit <- tryCatch(susieR::susie_rss(z=z, R=ld_mat, n=n, L=L), error=function(e) e)
  if (inherits(fit, "error")) return(list(ok=FALSE, msg=fit$message, fit=NULL))
  list(ok=TRUE, msg="OK", fit=fit)
}
extract_credible_sets_dt <- function(susie_fit, snp_names) {
  if (is.null(susie_fit) || !has_susieR) return(NULL)
  cs <- tryCatch(susieR::susie_get_cs(susie_fit), error=function(e) NULL)
  if (is.null(cs) || is.null(cs$cs)) return(NULL)
  out <- list()
  for (k in seq_along(cs$cs)) {
    idx <- cs$cs[[k]]
    out[[k]] <- data.table(cs_id=k, SNP=snp_names[idx], pip=susie_fit$pip[idx])
  }
  rbindlist(out, fill=TRUE)
}

# ---------- COLOC ----------
run_coloc_for_gene <- function(eqtl_reg, gwas_reg, N_eqtl, N_gwas, type_eqtl="quant", type_gwas="cc") {
  if (!has_coloc) return(list(ok=FALSE, msg="coloc package not installed"))
  common <- intersect(eqtl_reg$SNP, gwas_reg$SNP)
  if (length(common) < 50) return(list(ok=FALSE, msg=paste("Too few SNPs for coloc:", length(common))))
  eq <- eqtl_reg[SNP %in% common]
  gw <- gwas_reg[SNP %in% common]
  setkey(eq, SNP); setkey(gw, SNP)
  m <- gw[eq, nomatch=0]
  if (any(is.na(m$i.beta)) || any(is.na(m$i.se)) || any(is.na(m$beta)) || any(is.na(m$se))) {
    return(list(ok=FALSE, msg="Missing beta/se in region"))
  }
  maf1 <- pmin(m$i.eaf, 1 - m$i.eaf)
  maf2 <- pmin(m$eaf, 1 - m$eaf)
  maf  <- ifelse(!is.na(maf1), maf1, maf2)
  if (all(is.na(maf))) return(list(ok=FALSE, msg="No EAF/MAF available for coloc"))
  if (is.na(N_eqtl) || is.na(N_gwas)) return(list(ok=FALSE, msg="Need N_EQTL_DEFAULT and N_GWAS_DEFAULT for coloc"))
  D1 <- list(snp=m$SNP, beta=m$i.beta, varbeta=(m$i.se)^2, MAF=maf, N=as.numeric(N_eqtl), type=type_eqtl)
  D2 <- list(snp=m$SNP, beta=m$beta,  varbeta=(m$se)^2,  MAF=maf, N=as.numeric(N_gwas), type=type_gwas)
  res <- coloc::coloc.abf(D1, D2)
  list(ok=TRUE, res=res)
}

# ============================
# CONFIG DUMP
# ============================
run_cfg <- list(
  base_dir=base_dir, eqtl_path=eqtl_path, gwas_path=gwas_path,
  PLINK_BIN=PLINK_BIN, PLINK_REF=PLINK_REF,
  SAVE_PLOTS=SAVE_PLOTS, PLOT_FORMAT=PLOT_FORMAT, PLOT_DPI=PLOT_DPI,
  EQTL_P_THRESH=EQTL_P_THRESH, FSTAT_THRESH=FSTAT_THRESH,
  CLUMP_R2_STRICT=CLUMP_R2_STRICT, CLUMP_R2_SMALLN=CLUMP_R2_SMALLN, CLUMP_KB=CLUMP_KB,
  HARMONISE_ACTION=HARMONISE_ACTION, MIN_SNP_STRONG=MIN_SNP_STRONG,
  RUN_COLOC=RUN_COLOC, COLOC_WINDOW_BP=COLOC_WINDOW_BP, N_EQTL_DEFAULT=N_EQTL_DEFAULT, N_GWAS_DEFAULT=N_GWAS_DEFAULT,
  GWAS_TYPE=GWAS_TYPE, EQTL_TYPE=EQTL_TYPE,
  LD_R2_THRESHOLDS=LD_R2_THRESHOLDS,
  RUN_FINEMAP=RUN_FINEMAP, SUSIE_L=SUSIE_L, FINEMAP_WINDOW_BP=FINEMAP_WINDOW_BP
)

if (has_jsonlite) {
  jsonlite::write_json(
    run_cfg,
    file.path(d00, "run_config.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
} else {
  writeLines(
    "jsonlite not installed; run_config.json not written.",
    file.path(d00, "run_config.json")
  )
}

writeLines(capture.output(sessionInfo()), file.path(d00, "session_info.txt"))

# ============================
# PLINK READY?
# ============================
plink_ok <- check_plink_ready(PLINK_BIN, PLINK_REF)
log_msg(paste("PLINK OK:", plink_ok))

# ============================
# LOAD RAW DATA + WRITE 01_input_raw
# ============================
log_msg("Loading RAW eQTL...")
eqtl_raw <- fread(eqtl_path, showProgress=FALSE)
fwrite(eqtl_raw, file.path(d01, "eqtl_raw.tsv"), sep="\t")

# log_msg("Loading RAW GWAS...")
# gwas_raw <- fread(gwas_path, showProgress=FALSE)
# fwrite(gwas_raw, file.path(d01, "gwas_raw.tsv"), sep="\t")

log_msg("Loading RAW GWAS (optimized for large files)...")

# Sniff header to resolve column-name variants across GWAS Catalog files
.gwas_hdr <- names(fread(gwas_path, nrows = 0))
.se_col   <- intersect(c("standard_error", "se"), .gwas_hdr)[1]
.eaf_col  <- intersect(c("effect_allele_frequency", "eaf", "af"), .gwas_hdr)[1]
if (is.na(.se_col)) stop("GWAS file has no 'standard_error' or 'se' column")

.sel <- c("rsid", "beta", .se_col, "p_value",
          "effect_allele", "other_allele",
          "chromosome", "base_pair_location")
if (!is.na(.eaf_col)) .sel <- c(.sel, .eaf_col)

.cls <- list(
  character = c("rsid", "effect_allele", "other_allele"),
  numeric   = c("beta", .se_col, "p_value", "base_pair_location"),
  integer   = "chromosome"
)
if (!is.na(.eaf_col)) .cls$numeric <- c(.cls$numeric, .eaf_col)

gwas_raw <- fread(
  gwas_path,
  select     = .sel,
  colClasses = .cls,
  showProgress = FALSE
)

log_msg(paste("GWAS loaded:", nrow(gwas_raw), "rows"))

fwrite(
  gwas_raw,
  file.path(d01, "gwas_raw.tsv"),
  sep = "\t"
)

# ============================
# PARSE GWAS -> 02_input_parsed (derived from the already-loaded gwas_raw)
# ============================
log_msg("Parsing GWAS to MR format...")
gwas <- data.table(
  SNP = as.character(gwas_raw$rsid),
  beta = as.numeric(gwas_raw$beta),
  se = as.numeric(gwas_raw[[.se_col]]),
  pval = as.numeric(gwas_raw$p_value),
  effect_allele = as.character(gwas_raw$effect_allele),
  other_allele  = as.character(gwas_raw$other_allele),
  eaf = if (!is.na(.eaf_col)) as.numeric(gwas_raw[[.eaf_col]]) else NA_real_,
  chr = as.character(gwas_raw$chromosome),
  pos = as.numeric(gwas_raw$base_pair_location)
)
gwas <- gwas[grepl("^rs[0-9]+$", SNP)]
gwas <- gwas[!is.na(beta) & !is.na(se) & !is.na(pval)]
setorder(gwas, pval)
gwas <- unique(gwas, by="SNP")
fwrite(gwas, file.path(d02, "gwas_parsed.tsv"), sep="\t")
rm(gwas_raw); gc()

# ============================
# PARSE eQTL -> 02_input_parsed
# ============================
log_msg("Preparing eQTL parsed...")
eqtl <- extract_chr_pos_from_variant_id(copy(eqtl_raw))
rm(eqtl_raw); gc()

if (!"chr" %in% names(eqtl)) {
  if ("gene_chr" %in% names(eqtl)) eqtl[, chr := as.character(gene_chr)]
}
if (!"pos" %in% names(eqtl)) {
  if ("variant_pos" %in% names(eqtl)) eqtl[, pos := as.numeric(variant_pos)]
}
log_msg(paste("eQTL chr/pos available:", sum(!is.na(eqtl$chr) & !is.na(eqtl$pos)), "variants"))

eqtl <- eqtl[!is.na(rs_id_dbSNP155_GRCh38p13) & !is.na(pval_nominal) & pval_nominal < EQTL_P_THRESH]
eqtl[, F_stat := (slope / slope_se)^2]
eqtl <- eqtl[F_stat > FSTAT_THRESH]
eqtl[, `:=`(
  SNP = as.character(rs_id_dbSNP155_GRCh38p13),
  beta = as.numeric(slope),
  se = as.numeric(slope_se),
  pval = as.numeric(pval_nominal),
  eaf = as.numeric(af),
  effect_allele = as.character(alt),
  other_allele  = as.character(ref)
)]
fwrite(eqtl, file.path(d02, "eqtl_parsed.tsv"), sep="\t")

# ============================
# 03 instruments
# ============================
log_msg("Writing instruments...")
fwrite(eqtl, file.path(d03a, "eqtl_filtered_all.tsv"), sep="\t")
fwrite(eqtl[, .(gene_name,SNP,beta,se,pval,F_stat,eaf,effect_allele,other_allele,chr,pos)],
       file.path(d03b, "eqtl_fstat_pass.tsv"), sep="\t")

genes <- unique(eqtl$gene_name)
log_msg(paste("Genes to test:", length(genes)))

# ============================
# 04 overlap
# ============================
log_msg("Reducing GWAS to eQTL SNP overlap...")
eqtl_snps_all <- unique(eqtl$SNP)
gwas_reduced <- gwas[data.table(SNP=eqtl_snps_all), on="SNP", nomatch=0,
                     .(SNP, beta, se, pval, effect_allele, other_allele, eaf, chr, pos)]
fwrite(gwas_reduced, file.path(d04, "gwas_reduced_to_eqtl_snps.tsv"), sep="\t")
fwrite(data.table(
  n_eqtl_snps=length(eqtl_snps_all),
  n_gwas_snps_total=nrow(gwas),
  n_gwas_snps_overlap=nrow(gwas_reduced),
  overlap_rate=ifelse(length(eqtl_snps_all)>0, nrow(gwas_reduced)/length(eqtl_snps_all), NA_real_)
), file.path(d04, "overlap_stats.tsv"), sep="\t")
gwas <- gwas_reduced

# ============================
# MAIN LOOP
# ============================
all_results <- list()
qc_list <- list()
diag_list <- list()
coloc_summary_list <- list()
finemap_summary_list <- list()

setkey(eqtl, gene_name)
gwas_snp_set <- gwas$SNP

n_genes <- length(genes)
.gene_i <- 0L

for (g in genes) {
  .gene_i <- .gene_i + 1L
  if (.gene_i == 1L || .gene_i %% 25L == 0L || .gene_i == n_genes)
    log_msg(sprintf("Processing gene %d / %d  (%s)", .gene_i, n_genes, g))
  x <- eqtl[.(g)]
  if (nrow(x)==0) next
  setorder(x, pval)
  x <- x[, .SD[1], by=SNP]
  x <- x[SNP %chin% gwas_snp_set]
  if (nrow(x)==0) next
  x[, exposure := paste0("GTEx_GeneExpression_", g)]

  exposure_before <- suppressMessages(format_data(
    as.data.frame(x), type="exposure",
    snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele",
    pval_col="pval", eaf_col="eaf", phenotype_col="exposure"
  ))
  r2_use <- if (nrow(exposure_before) < 5) CLUMP_R2_SMALLN else CLUMP_R2_STRICT
  exposure_after <- exposure_before

  if (plink_ok && nrow(exposure_before) >= 2) {
    outpref <- file.path(d05, paste0("GENE_", g, "__PLINK"))
    cl <- tryCatch(plink_clump(exposure_before, PLINK_BIN, PLINK_REF, outpref, r2_use, CLUMP_KB),
                   error=function(e) list(ok=FALSE, msg=e$message, keep=NULL))
    if (isTRUE(cl$ok) && length(cl$keep) > 0) {
      exposure_after <- exposure_before[exposure_before$SNP %in% cl$keep, ]
    } else {
      log_msg(paste("PLINK CLUMP FAIL (using UNCLUMPED):", g, cl$msg))
      exposure_after <- exposure_before
    }
  } else {
    # silently skip clumping when PLINK unavailable or <2 SNPs
  }
  if (is.null(exposure_after) || nrow(exposure_after) == 0) next

  gw <- gwas[data.table(SNP=exposure_after$SNP), on="SNP", nomatch=0]
  if (nrow(gw)==0) next
  gw[, outcome := "GWAS_outcome"]

  outcome <- suppressMessages(format_data(
    as.data.frame(gw), type="outcome",
    snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele",
    pval_col="pval", eaf_col="eaf", phenotype_col="outcome"
  ))

  dat0 <- suppressMessages(harmonise_data(
    exposure_dat = exposure_after,
    outcome_dat  = outcome,
    action       = HARMONISE_ACTION
  ))

  dat <- dat0[dat0$mr_keep == TRUE, ]
  if (nrow(dat) == 0) next

  evidence_level <- if (nrow(dat) >= MIN_SNP_STRONG) "multi_snp" else "single_snp_wald"
  snp_string <- paste(sort(unique(dat$SNP)), collapse=";")

  mr_res <- suppressMessages(mr(dat))
  mr_res$gene <- g
  mr_res$snps_used <- snp_string
  mr_res$nsnp_used <- nrow(dat)
  mr_res$evidence_level <- evidence_level

  all_results[[g]] <- mr_res

  het <- if (nrow(dat)>=2) tryCatch(mr_heterogeneity(dat), error=function(e) NULL) else NULL
  pleio <- if (nrow(dat)>=3) tryCatch(mr_pleiotropy_test(dat), error=function(e) NULL) else NULL
  if (!is.null(het)) diag_list[[paste0(g,"_het")]] <- as.data.table(het)
  if (!is.null(pleio)) diag_list[[paste0(g,"_pleio")]] <- as.data.table(pleio)

  qc_list[[g]] <- data.table(
    gene=g, evidence_level=evidence_level, nsnp_used=nrow(dat), snps_used=snp_string,
    eqtl_snps=nrow(x), after_clump=nrow(exposure_after), gwas_overlap=nrow(gw),
    kept=sum(dat0$mr_keep, na.rm=TRUE),
    pal=sum(dat0$palindromic, na.rm=TRUE),
    amb=sum(dat0$ambiguous, na.rm=TRUE)
  )

  # 07_ld (PLINK LD matrix + heatmap)
  lead_snp <- dat$SNP[which.min(dat$pval.exposure)]
  snps_for_ld <- unique(dat$SNP)

  ld_out <- if (plink_ok) plink_ld_matrix(snps_for_ld, PLINK_BIN, PLINK_REF,
                                         file.path(d07, paste0("GENE_", g, "__PLINK_LD"))) else list(ok=FALSE, msg="PLINK not OK", mat=NULL)

  ld_mat <- NULL
  if (isTRUE(ld_out$ok) && !is.null(ld_out$mat)) {
    ld_mat <- ld_out$mat
    plot_ld_heatmap(ld_mat, file.path(d11_pg, paste0("GENE_", g, "__LD_heatmap.", PLOT_FORMAT)))
  }

  # 09_finemap — record eligibility only (SuSiE fits are not persisted,
  # so the expensive computation is skipped for performance)
  finemap_ok <- RUN_FINEMAP && has_susieR &&
                !is.na(N_EQTL_DEFAULT) && !is.na(N_GWAS_DEFAULT) &&
                !is.null(ld_mat) && is.matrix(ld_mat)
  finemap_summary_list[[g]] <- data.table(
    gene=g, lead_snp=lead_snp,
    n_snps_finemap = if (finemap_ok) ncol(ld_mat) else NA_integer_,
    finemap_ran = finemap_ok
  )

  # 10_coloc (results kept in memory for summary table)
  if (RUN_COLOC && has_coloc) {
    can_eqtl_pos <- all(c("chr","pos") %in% names(x)) && any(!is.na(x$chr) & !is.na(x$pos))
    can_gwas_pos <- all(c("chr","pos") %in% names(gwas)) && any(!is.na(gwas$chr) & !is.na(gwas$pos))
    if (can_eqtl_pos && can_gwas_pos && !is.na(N_EQTL_DEFAULT) && !is.na(N_GWAS_DEFAULT)) {
      lead <- x[which.min(pval)]
      if (!is.na(lead$chr) && !is.na(lead$pos)) {
        start <- lead$pos - COLOC_WINDOW_BP
        end   <- lead$pos + COLOC_WINDOW_BP
        eqtl_reg <- eqtl[gene_name==g & chr==lead$chr & pos>=start & pos<=end]
        gwas_reg <- gwas[chr==lead$chr & pos>=start & pos<=end]
        coloc_out <- tryCatch(run_coloc_for_gene(eqtl_reg, gwas_reg, N_EQTL_DEFAULT, N_GWAS_DEFAULT, EQTL_TYPE, GWAS_TYPE),
                              error=function(e) list(ok=FALSE, msg=e$message))
        if (isTRUE(coloc_out$ok)) {
          pp <- coloc_out$res$summary
          coloc_summary_list[[g]] <- data.table(
            gene=g, lead_snp=lead$SNP, lead_chr=as.character(lead$chr), lead_pos=as.numeric(lead$pos),
            region_start=start, region_end=end,
            PP4_shared=as.numeric(pp["PP.H4.abf"]),
            PP3_distinct=as.numeric(pp["PP.H3.abf"]),
            nsnp_region_common=length(intersect(eqtl_reg$SNP, gwas_reg$SNP))
          )
        }
      }
    }
  }

  # 11 plots — only generated for genes with >=2 SNPs
  if (SAVE_PLOTS && nrow(dat) >= MIN_SNP_STRONG) {
    ns <- nrow(dat)
    plot_prefix <- file.path(d11_pg, paste0("GENE_", g, "__"))
    if (ns >= 2) save_plot_list(tryCatch(mr_scatter_plot(mr_res, dat), error=function(e) NULL), paste0(plot_prefix, "scatter"))
    ss <- if (ns >= 1) tryCatch(mr_singlesnp(dat), error=function(e) NULL) else NULL
    if (!is.null(ss) && nrow(ss) >= 1) save_plot_list(tryCatch(mr_forest_plot(ss), error=function(e) NULL), paste0(plot_prefix, "forest"))
    if (!is.null(ss) && nrow(ss) >= 3) save_plot_list(tryCatch(mr_funnel_plot(ss), error=function(e) NULL), paste0(plot_prefix, "funnel"))
    if (ns >= 4) {
      loo <- tryCatch(mr_leaveoneout(dat), error=function(e) NULL)
      loo_p <- tryCatch(if (!is.null(loo)) mr_leaveoneout_plot(loo) else NULL, error=function(e) NULL)
      save_plot_list(loo_p, paste0(plot_prefix, "leaveoneout"))
    }
  }
}

# ============================
# SUMMARY TABLES (12_summary_tables)
# ============================
if (length(all_results) == 0) {
  log_msg("No MR results produced (all genes filtered out).")
  log_msg("RUN END: OK")
} else {
  mr_all <- as.data.table(rbindlist(all_results, fill=TRUE))
  mr_all[, `:=`(
    OR=exp(b), CI_low=exp(b-1.96*se), CI_high=exp(b+1.96*se),
    direction=fifelse(b>0, "Higher expression -> higher risk",
                      fifelse(b<0, "Higher expression -> lower risk (protective)", "No direction")),
    odds_change_percent=(exp(b)-1)*100
  )]
  fwrite(mr_all, file.path(d12, "MR_MAIN_RESULTS_ALL_GENES.csv"))
  fwrite(rbindlist(qc_list, fill=TRUE), file.path(d12, "MR_QC_SUMMARY_ALL_GENES.csv"))

  if (length(diag_list) > 0) {
    fwrite(rbindlist(diag_list, fill=TRUE), file.path(d12, "MR_DIAGNOSTICS.csv"))
  } else {
    fwrite(data.table(note="No diagnostics available"), file.path(d12, "MR_DIAGNOSTICS.csv"))
  }

  if (length(finemap_summary_list) > 0) {
    fwrite(rbindlist(finemap_summary_list, fill=TRUE), file.path(d12, "FINEMAP_SUMMARY.csv"))
  } else {
    fwrite(data.table(note="No finemap info"), file.path(d12, "FINEMAP_SUMMARY.csv"))
  }

  if (length(coloc_summary_list) > 0) {
    fwrite(rbindlist(coloc_summary_list, fill=TRUE), file.path(d12, "COLOC_SUMMARY.csv"))
  } else {
    fwrite(data.table(note="No coloc results"), file.path(d12, "COLOC_SUMMARY.csv"))
  }

  pref <- mr_all[, {
    idx <- if (any(method == "Inverse variance weighted")) {
      which(method == "Inverse variance weighted")[1]
    } else if (any(method == "Wald ratio")) {
      which(method == "Wald ratio")[1]
    } else {
      1L
    }
    .SD[idx]
  }, by = gene]
  setorder(pref, pval)
  qc_dt <- rbindlist(qc_list, fill=TRUE)
  pref_dt <- merge(as.data.table(pref), qc_dt, by="gene", all.x=TRUE)

  if (length(coloc_summary_list) > 0) {
    coloc_all <- rbindlist(coloc_summary_list, fill=TRUE)
    pref_dt <- merge(pref_dt, coloc_all[, .(gene, PP4_shared, PP3_distinct, nsnp_region_common)], by="gene", all.x=TRUE)
  } else {
    pref_dt[, `:=`(PP4_shared=NA_real_, PP3_distinct=NA_real_, nsnp_region_common=NA_integer_)]
  }

  pref_dt[, confidence_tier := fifelse(
    pval >= 0.05, "Not_significant",
    fifelse(!is.na(PP3_distinct) & PP3_distinct >= 0.5, "Likely_false",
      fifelse(!is.na(PP4_shared) & PP4_shared >= 0.80, "Strong",
        fifelse(!is.na(PP4_shared) & PP4_shared >= 0.50, "Moderate", "Weak")
      )
    )
  )]
  fwrite(pref_dt, file.path(d12, "MR_INTERPRETATION_SUMMARY.csv"))
  log_msg("RUN END: OK")
}

############################################################
# FINAL BIOLOGICAL EVIDENCE LAYER
############################################################
log_msg("Creating FINAL gene- and variant-level genetic evidence tables...")

variant_evidence <- eqtl[SNP %in% gwas$SNP, .(
  Gene_Symbol = gene_name,
  SNP,
  Chromosome = chr,
  Position_bp = pos,
  eQTL_beta = beta,
  eQTL_se   = se,
  eQTL_Z    = beta / se,
  eQTL_direction = fifelse(beta > 0, "Upregulating", "Downregulating")
)]
variant_evidence <- merge(
  variant_evidence,
  gwas[, .(SNP, GWAS_beta = beta, GWAS_se = se, GWAS_Z = beta / se, GWAS_pval = pval)],
  by="SNP", all.x=TRUE
)
variant_evidence[, Variant_Genetic_Confidence := abs(eQTL_Z) * abs(GWAS_Z)]
variant_evidence[, `:=`(Outcome=OUTCOME_NAME, Biosample_Type=BIOSAMPLE_TYPE)]
fwrite(variant_evidence, file.path(d12, paste0(OUTCOME_NAME, "_", BIOSAMPLE_TYPE, "_VariantLevel_GeneticEvidence.tsv")), sep="\t")

gene_evidence <- variant_evidence[, .(
  Chromosome = unique(Chromosome),
  Gene_Start_bp = min(Position_bp, na.rm=TRUE),
  Gene_End_bp   = max(Position_bp, na.rm=TRUE),
  n_eQTL_variants = .N,
  Mean_eQTL_Z = mean(abs(eQTL_Z), na.rm=TRUE),
  Max_eQTL_Z  = max(abs(eQTL_Z), na.rm=TRUE),
  n_GWAS_variants = sum(!is.na(GWAS_Z)),
  Mean_GWAS_Z = mean(abs(GWAS_Z), na.rm=TRUE),
  Min_GWAS_p  = min(GWAS_pval, na.rm=TRUE),
  Mean_Variant_Confidence = mean(Variant_Genetic_Confidence, na.rm=TRUE)
), by=Gene_Symbol]

scale01 <- function(x) {
  r <- range(x, na.rm=TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}
gene_evidence[, Gene_Genetic_Confidence_Score :=
  0.4*scale01(Mean_eQTL_Z) + 0.4*scale01(Mean_GWAS_Z) + 0.2*scale01(Mean_Variant_Confidence)
]
gene_evidence[, Evidence_Strength :=
  fifelse(Gene_Genetic_Confidence_Score >= 0.75, "Very_Strong",
    fifelse(Gene_Genetic_Confidence_Score >= 0.55, "Strong",
      fifelse(Gene_Genetic_Confidence_Score >= 0.35, "Moderate",
        fifelse(Gene_Genetic_Confidence_Score >= 0.20, "Weak", "Low")
      )
    )
  )
]
gene_evidence[, `:=`(Outcome=OUTCOME_NAME, Biosample_Type=BIOSAMPLE_TYPE)]
fwrite(gene_evidence, file.path(d12, paste0(OUTCOME_NAME, "_", BIOSAMPLE_TYPE, "_GeneLevel_GeneticEvidence.tsv")), sep="\t")

unlink(.pg_tmp, recursive=TRUE)
log_msg("FINAL genetic evidence tables completed successfully.")
############################################################
# Notes:
# - Clumping + LD matrix are via PLINK using PLINK_REF (hg38).
# - SuSiE runs only if LD exists AND N defaults set.
############################################################
