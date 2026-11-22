# 07_sensitivity_outliers.R
# Goal: finish the outlier follow-up after the QC step.
#   Part A: repeat the DESeq2 analysis without each flagged patient (P3, P17, P18) and
#           without all three together, and compare each result with the main analysis.
#           Optional: leave-one-patient-out over all 18 patients (slow, off by default).
#   Part B: check whether the "normal" sample of P17 (and of P3, P18) looks like tumour
#           tissue, using marker genes (tumour-up genes and normal-epithelium genes).
# Input : data/processed/counts_filtered.csv, data/processed/vst_matrix.csv,
#         metadata/sample_metadata_clean.csv, metadata/gene_annotation.csv,
#         results/deseq2_all_genes.csv (main analysis from script 04)
# Output: results/outlier_sensitivity_summary.csv, results/outlier_sensitivity_summary.md,
#         results/deseq2_sensitivity_without_<patients>.csv,
#         results/marker_tumour_likeness.csv, results/marker_genes_flagged_normals.csv,
#         figures/outlier_sensitivity_scatter.png, figures/marker_tumour_likeness.png

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
set.seed(1)

# ---- settings for this script ----
FLAGGED_PATIENTS  <- c("P3", "P17", "P18")  # patients flagged in QC (see results/qc_outlier_flags.csv)
RUN_LEAVE_ONE_OUT <- FALSE                  # TRUE = also remove each of the 18 patients in turn (slow)

up_lab   <- paste("Up in", GROUP_TEST)
down_lab <- paste("Down in", GROUP_TEST)

# ------------------------------------------------------------------ 1. inputs
meta <- read.csv("metadata/sample_metadata_clean.csv", stringsAsFactors = FALSE)
meta$group   <- factor(meta$group, levels = c(GROUP_REF, GROUP_TEST))
meta$patient <- factor(meta$patient)
rownames(meta) <- meta$sample_id

mat <- as.matrix(read.csv("data/processed/counts_filtered.csv", row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(mat), meta$sample_id))
storage.mode(mat) <- "integer"

vst <- as.matrix(read.csv("data/processed/vst_matrix.csv", row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(vst), meta$sample_id))

ann <- read.csv("metadata/gene_annotation.csv", colClasses = "character")
if (!"gene_name" %in% names(ann)) ann$gene_name <- NA_character_
ann <- ann[!duplicated(ann$gene_id), ]

main <- read.csv("results/deseq2_all_genes.csv", colClasses = c(gene_id = "character"))
stopifnot(all(c("gene_id", "log2FC_shrunk", "direction") %in% names(main)))

missing_flagged <- setdiff(FLAGGED_PATIENTS, as.character(meta$patient))
if (length(missing_flagged) > 0)
  stop("Flagged patients not found in the metadata: ", paste(missing_flagged, collapse = ", "))

# --------------------------------------------------------------- 2. functions
# Same model as script 04: design ~ patient + group, apeglm shrinkage, same thresholds.
run_deseq <- function(counts, md) {
  md$group   <- factor(md$group, levels = c(GROUP_REF, GROUP_TEST))
  md$patient <- droplevels(factor(md$patient))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = md, design = ~ patient + group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", GROUP_TEST, GROUP_REF), alpha = PADJ_CUTOFF)
  coef_name <- paste0("group_", GROUP_TEST, "_vs_", GROUP_REF)
  stopifnot(coef_name %in% resultsNames(dds))
  shr <- tryCatch(
    lfcShrink(dds, coef = coef_name, res = res, type = "apeglm"),
    error = function(e) {
      message("apeglm shrinkage unavailable (", conditionMessage(e), "): unshrunk values used.")
      NULL
    })
  list(dds = dds, res = res, shr = shr)
}

make_table <- function(fit) {
  res <- fit$res
  tab <- data.frame(
    gene_id       = rownames(res),
    baseMean      = res$baseMean,
    log2FC        = res$log2FoldChange,
    log2FC_shrunk = if (!is.null(fit$shr)) fit$shr$log2FoldChange else res$log2FoldChange,
    stat          = res$stat,
    pvalue        = res$pvalue,
    padj          = res$padj,
    stringsAsFactors = FALSE
  )
  i <- match(tab$gene_id, ann$gene_id)
  tab$symbol <- ann$symbol[i]
  tab$label  <- ifelse(is.na(tab$symbol) | tab$symbol == "", tab$gene_id, tab$symbol)
  sig <- !is.na(tab$padj) & tab$padj < PADJ_CUTOFF & abs(tab$log2FC_shrunk) >= LFC_CUTOFF
  tab$direction <- ifelse(sig & tab$log2FC_shrunk > 0, up_lab,
                   ifelse(sig & tab$log2FC_shrunk < 0, down_lab, "Not significant"))
  tab <- tab[order(tab$padj, na.last = TRUE), ]
  rownames(tab) <- NULL
  tab
}

run_without <- function(patients) {
  keep <- !(as.character(meta$patient) %in% patients)
  fit_s <- run_deseq(mat[, keep, drop = FALSE], meta[colnames(mat)[keep], ])
  list(tab = make_table(fit_s), n_samples = sum(keep))
}

compare_to_main <- function(tab_s, scenario, excluded, n_samples) {
  m <- merge(main[, c("gene_id", "log2FC_shrunk", "direction")],
             tab_s[, c("gene_id", "log2FC_shrunk", "direction")],
             by = "gene_id", suffixes = c("_all", "_sens"))
  r      <- cor(m$log2FC_shrunk_all, m$log2FC_shrunk_sens, use = "complete.obs")
  s_all  <- m$gene_id[m$direction_all  != "Not significant"]
  s_sens <- m$gene_id[m$direction_sens != "Not significant"]
  both   <- intersect(s_all, s_sens)
  uni    <- union(s_all, s_sens)
  data.frame(
    scenario = scenario, excluded = excluded, samples = n_samples,
    correlation = round(r, 3),
    significant_main = length(s_all), significant_scenario = length(s_sens),
    shared = length(both),
    pct_of_main_retained = round(100 * length(both) / length(s_all), 1),
    pct_of_union = round(100 * length(both) / length(uni), 1),
    stringsAsFactors = FALSE)
}

# ------------------------------------------------ 3. Part A: sensitivity analyses
scenarios <- list()
for (p in FLAGGED_PATIENTS) scenarios[[paste("without", p)]] <- p
scenarios[[paste("without", paste(FLAGGED_PATIENTS, collapse = ", "))]] <- FLAGGED_PATIENTS

summary_rows <- list()
scatter_data <- list()
for (nm in names(scenarios)) {
  pat <- scenarios[[nm]]
  message("Running scenario: ", nm, " ...")
  out <- run_without(pat)
  write.csv(out$tab, paste0("results/deseq2_sensitivity_without_", paste(pat, collapse = "_"), ".csv"),
            row.names = FALSE)
  summary_rows[[nm]] <- compare_to_main(out$tab, nm, paste(pat, collapse = ", "), out$n_samples)
  m <- merge(main[, c("gene_id", "log2FC_shrunk")], out$tab[, c("gene_id", "log2FC_shrunk")],
             by = "gene_id", suffixes = c("_all", "_sens"))
  m$scenario <- nm
  scatter_data[[nm]] <- m
}

if (RUN_LEAVE_ONE_OUT) {
  for (p in levels(meta$patient)) {
    message("Leave-one-patient-out: ", p, " ...")
    out <- run_without(p)
    summary_rows[[paste("LOO", p)]] <- compare_to_main(out$tab, paste("leave out", p), p, out$n_samples)
  }
}

summ <- do.call(rbind, summary_rows)
rownames(summ) <- NULL
write.csv(summ, "results/outlier_sensitivity_summary.csv", row.names = FALSE)

sd_all <- do.call(rbind, scatter_data)
sd_all$scenario <- factor(sd_all$scenario, levels = names(scenarios))
r_lab <- summ[match(levels(sd_all$scenario), summ$scenario), ]
facet_names <- setNames(paste0(r_lab$scenario, " (r = ", sprintf("%.3f", r_lab$correlation), ")"),
                        r_lab$scenario)
p_s <- ggplot(sd_all, aes(log2FC_shrunk_all, log2FC_shrunk_sens)) +
  geom_point(alpha = 0.25, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, colour = "#D55E00", linetype = "dashed") +
  facet_wrap(~ scenario, labeller = labeller(scenario = facet_names)) +
  labs(x = "log2FC, all samples", y = "log2FC, scenario",
       title = "Sensitivity to the QC-flagged patients") +
  theme_minimal(base_size = 11)
ggsave("figures/outlier_sensitivity_scatter.png", p_s, width = 8, height = 7, dpi = 200)

md_tab <- c("| Scenario | Samples | r (log2FC) | Significant (main -> scenario) | Shared | Main genes retained | Overlap of union |",
            "|---|---|---|---|---|---|---|",
            paste0("| ", summ$scenario, " | ", summ$samples, " | ", sprintf("%.3f", summ$correlation),
                   " | ", summ$significant_main, " -> ", summ$significant_scenario,
                   " | ", summ$shared, " | ", summ$pct_of_main_retained, "% | ", summ$pct_of_union, "% |"))

# ------------------------------------------------ 4. Part B: marker genes in normal samples
markers <- data.frame(
  class  = c(rep("Tumour-up genes", 8), rep("Normal-epithelium genes", 6), rep("Proliferation genes", 2)),
  symbol = c("MMP7", "KRT80", "COL11A1", "CST1", "INHBA", "ETV4", "CLDN1", "CDH3",
             "BEST4", "CA7", "OTOP2", "OTOP3", "UNC5C", "ABCG2",
             "MKI67", "TOP2A"),
  stringsAsFactors = FALSE)
markers$gene_id <- ann$gene_id[match(markers$symbol, ann$symbol)]
absent <- markers$symbol[is.na(markers$gene_id) | !(markers$gene_id %in% rownames(vst))]
if (length(absent) > 0) message("Marker genes not available (skipped): ", paste(absent, collapse = ", "))
markers <- markers[!is.na(markers$gene_id) & markers$gene_id %in% rownames(vst), ]
stopifnot(nrow(markers) >= 4)

normal_ids <- meta$sample_id[meta$group == GROUP_REF]
tumor_ids  <- meta$sample_id[meta$group == GROUP_TEST]

# "Tumour-likeness" of a sample for each marker gene:
#   0 = same level as the other normal samples, 1 = same level as the tumour average.
# The sample itself is left out of the normal average, so a normal sample is not compared with itself.
likeness <- function(s) {
  others <- setdiff(normal_ids, s)
  mu_n <- rowMeans(vst[markers$gene_id, others, drop = FALSE])
  mu_t <- rowMeans(vst[markers$gene_id, tumor_ids, drop = FALSE])
  (vst[markers$gene_id, s] - mu_n) / (mu_t - mu_n)
}
lk <- sapply(meta$sample_id, likeness)
rownames(lk) <- markers$symbol

main_markers <- markers$class != "Proliferation genes"
score <- apply(lk[main_markers, , drop = FALSE], 2, median)
sc <- data.frame(sample_id = names(score), patient = as.character(meta[names(score), "patient"]),
                 group = as.character(meta[names(score), "group"]),
                 tumour_likeness = round(unname(score), 3), stringsAsFactors = FALSE)
sc$flagged <- sc$patient %in% FLAGGED_PATIENTS
sc <- sc[order(sc$group, -sc$tumour_likeness), ]
rownames(sc) <- NULL
write.csv(sc, "results/marker_tumour_likeness.csv", row.names = FALSE)

# gene-by-gene detail for the flagged normal samples
detail <- do.call(rbind, lapply(FLAGGED_PATIENTS, function(p) {
  s <- meta$sample_id[meta$patient == p & meta$group == GROUP_REF]
  others <- setdiff(normal_ids, s)
  data.frame(patient = p, sample_id = s, class = markers$class, symbol = markers$symbol,
             mean_other_normals = round(rowMeans(vst[markers$gene_id, others, drop = FALSE]), 2),
             mean_tumours       = round(rowMeans(vst[markers$gene_id, tumor_ids, drop = FALSE]), 2),
             this_normal_sample = round(vst[markers$gene_id, s], 2),
             tumour_likeness    = round(lk[markers$symbol, s], 2), stringsAsFactors = FALSE)
}))
write.csv(detail, "results/marker_genes_flagged_normals.csv", row.names = FALSE)

pd <- sc
pd$label <- ifelse(pd$flagged & pd$group == GROUP_REF, pd$patient, "")
pd$sample_id <- factor(pd$sample_id, levels = pd$sample_id[order(pd$tumour_likeness)])
p_m <- ggplot(pd, aes(tumour_likeness, sample_id, colour = group)) +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey50") +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = -0.4, size = 3.5, colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = setNames(c("#0072B2", "#D55E00"), c(GROUP_REF, GROUP_TEST))) +
  labs(x = "Tumour-likeness of marker genes (0 = typical normal, 1 = tumour average)",
       y = NULL, colour = NULL,
       title = "Do the QC-flagged normal samples look like tumour tissue?") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())
ggsave("figures/marker_tumour_likeness.png", p_m, width = 8, height = 6, dpi = 200)

nrm <- sc[sc$group == GROUP_REF, ]
others_med <- median(nrm$tumour_likeness[!nrm$flagged])
flag_txt <- paste0("- ", nrm$patient[nrm$flagged], " (normal): tumour-likeness ",
                   sprintf("%.2f", nrm$tumour_likeness[nrm$flagged]))

# ------------------------------------------------ 5. summary file
writeLines(c(
  "# Outlier follow-up (script 07)", "",
  paste0("Flagged patients: ", paste(FLAGGED_PATIENTS, collapse = ", "),
         ". Same model and thresholds as script 04 (design ~ patient + group, padj < ", PADJ_CUTOFF,
         ", |shrunken log2FC| >= ", LFC_CUTOFF, ")."), "",
  "## A. Sensitivity analyses", md_tab, "",
  "Interpretation: a high correlation and a large share of the main genes retained mean that the",
  "conclusions do not depend on these patients.", "",
  "## B. Marker genes in the normal samples", "",
  paste0("Median tumour-likeness of the other (unflagged) normal samples: ", sprintf("%.2f", others_med)),
  flag_txt, "",
  "How to read it: 0 = the sample behaves like the other normal samples, 1 = like the average tumour.",
  "A flagged normal sample with a clearly higher value than the others may contain tumour tissue;",
  "a value close to the others points to a technical cause (library quality, depth) instead.",
  "Markers used: tumour-up (MMP7, KRT80, COL11A1, CST1, INHBA, ETV4, CLDN1, CDH3) and",
  "normal-epithelium (BEST4, CA7, OTOP2, OTOP3, UNC5C, ABCG2). This is a check, not a proof."
), "results/outlier_followup_summary.md")

writeLines(capture.output(sessionInfo()), "metadata/sessionInfo_outliers.txt")
message(paste(md_tab, collapse = "\n"))
message("\nMarker tumour-likeness (normal samples):")
print(nrm)
