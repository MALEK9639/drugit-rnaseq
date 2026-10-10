# 00_config.R -- shared settings, sourced by every script.
# Run all scripts from the project root (the folder that contains README.md).

# ---- Dataset (change here only if the supervisor approves another dataset) ----
GSE        <- "GSE50760"
GENOME_TAG <- "GRCh38.p13"

# ---- Comparison (reference level first) ----
GROUP_REF  <- "normal"   # reference / control group
GROUP_TEST <- "tumor"    # group compared to the reference

# ---- Quality-control parameters (documented in the report) ----
MIN_COUNT  <- 10         # a gene is "expressed" in a sample if count >= MIN_COUNT
# Filtering rule: keep genes with count >= MIN_COUNT in at least
# (size of the smallest group) samples.
N_TOP_VAR  <- 500        # most variable genes used for the PCA

# ---- Folders ----
for (d in c("data/raw", "data/processed", "metadata", "figures", "results", "report"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists("scripts/00_config.R"))
  stop("Run scripts from the project root: setwd('path/to/drugit-rnaseq')")
