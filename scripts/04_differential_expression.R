# 04_differential_expression.R
# Goal: paired differential expression, tumour vs normal, with DESeq2
#       (design ~ patient + group), plus a sensitivity analysis without one patient.
# Input : data/processed/counts_filtered.csv, data/processed/vst_matrix.csv,
#         metadata/sample_metadata_clean.csv, metadata/gene_annotation.csv
# Output: results/deseq2_all_genes.csv, results/deseq2_significant_genes.csv,
#         results/deseq2_summary.md, results/deseq2_sensitivity_without_<patient>.csv,
#         results/sensitivity_summary.md, data/processed/normalized_counts.csv,
#         figures/pvalue_histogram.png, dispersion_plot.png, ma_plot.png, volcano_plot.png,
#         top_genes_heatmap.png, top_genes_paired.png, sensitivity_lfc_scatter.png

source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)
})
set.seed(1)

up_lab   <- paste("Up in", GROUP_TEST)
down_lab <- paste("Down in", GROUP_TEST)
cols <- setNames(c("#D55E00", "#0072B2", "grey70"), c(up_lab, down_lab, "Not significant"))
grp_cols <- setNames(c("#0072B2", "#D55E00"), c(GROUP_REF, GROUP_TEST))

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

if (file.exists("metadata/gene_annotation.csv")) {
  ann <- read.csv("metadata/gene_annotation.csv", colClasses = "character")
} else {
  ann <- data.frame(gene_id = rownames(mat), symbol = NA_character_, gene_name = NA_character_)
}
if (!"gene_name" %in% names(ann)) ann$gene_name <- NA_character_
ann <- ann[!duplicated(ann$gene_id), ]

message("Samples: ", ncol(mat), " | genes tested: ", nrow(mat),
        " | patients: ", nlevels(meta$patient))

# --------------------------------------------------------------- 2. functions
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
    lfcSE         = res$lfcSE,
    stat          = res$stat,
    pvalue        = res$pvalue,
    padj          = res$padj,
    stringsAsFactors = FALSE
  )
  i <- match(tab$gene_id, ann$gene_id)
  tab$symbol    <- ann$symbol[i]
  tab$gene_name <- ann$gene_name[i]
  tab$label     <- ifelse(is.na(tab$symbol) | tab$symbol == "", tab$gene_id, tab$symbol)
  sig <- !is.na(tab$padj) & tab$padj < PADJ_CUTOFF & abs(tab$log2FC_shrunk) >= LFC_CUTOFF
  tab$direction <- ifelse(sig & tab$log2FC_shrunk > 0, up_lab,
                   ifelse(sig & tab$log2FC_shrunk < 0, down_lab, "Not significant"))
  tab <- tab[order(tab$padj, na.last = TRUE), ]
  rownames(tab) <- NULL
  tab
}

fmt_genes <- function(d) {
  if (nrow(d) == 0) return("- none")
  paste0("- ", d$label, " (log2FC ", sprintf("%.2f", d$log2FC_shrunk),
         ", padj ", sprintf("%.1e", d$padj), ")")
}

# ------------------------------------------------------- 3. main analysis (n = 36)
fit <- run_deseq(mat, meta[colnames(mat), ])
tab <- make_table(fit)
sig_tab <- tab[tab$direction != "Not significant", ]

write.csv(tab, "results/deseq2_all_genes.csv", row.names = FALSE)
write.csv(sig_tab, "results/deseq2_significant_genes.csv", row.names = FALSE)
nc <- counts(fit$dds, normalized = TRUE)
write.csv(nc, "data/processed/normalized_counts.csv")

n_tested <- sum(!is.na(tab$padj))
n_up   <- sum(tab$direction == up_lab)
n_down <- sum(tab$direction == down_lab)
message("Tested: ", n_tested, " | ", up_lab, ": ", n_up, " | ", down_lab, ": ", n_down)

summ <- c(
  "# Differential expression summary", "",
  paste0("- Dataset: ", GSE, "; comparison: ", GROUP_TEST, " vs ", GROUP_REF, " (reference)"),
  paste0("- Design: ~ patient + group (paired); samples: ", ncol(mat), "; patients: ", nlevels(meta$patient)),
  paste0("- Genes tested (after filtering and DESeq2 independent filtering): ", n_tested, " of ", nrow(mat)),
  paste0("- Thresholds: padj < ", PADJ_CUTOFF, " and |shrunken log2FC| >= ", LFC_CUTOFF),
  paste0("- ", up_lab, ": ", n_up, " genes; ", down_lab, ": ", n_down, " genes"), "",
  paste0("## Top genes ", tolower(up_lab)),
  fmt_genes(head(tab[tab$direction == up_lab, ], 10)), "",
  paste0("## Top genes ", tolower(down_lab)),
  fmt_genes(head(tab[tab$direction == down_lab, ], 10))
)
writeLines(summ, "results/deseq2_summary.md")

# ------------------------------------------------------------------ 4. figures
# 4a. p-value histogram (diagnostic: a flat bottom plus a spike near 0 is healthy)
p_hist <- ggplot(tab[!is.na(tab$pvalue), ], aes(pvalue)) +
  geom_histogram(bins = 50, fill = "#4C78A8", colour = "white") +
  labs(x = "p-value", y = "Number of genes", title = "Distribution of raw p-values") +
  theme_minimal(base_size = 11)
ggsave("figures/pvalue_histogram.png", p_hist, width = 6, height = 4, dpi = 200)

# 4b. dispersion estimates (diagnostic of the model fit)
png("figures/dispersion_plot.png", width = 1400, height = 1100, res = 200)
plotDispEsts(fit$dds, main = "DESeq2 dispersion estimates")
invisible(dev.off())

# 4c. MA plot
p_ma <- ggplot(tab[!is.na(tab$log2FC_shrunk), ],
               aes(log10(baseMean + 1), log2FC_shrunk, colour = direction)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = cols) +
  labs(x = "log10 mean normalised count + 1", y = "log2 fold change (shrunken)",
       colour = NULL, title = paste0(GSE, ": MA plot (", GROUP_TEST, " vs ", GROUP_REF, ")")) +
  theme_minimal(base_size = 11)
ggsave("figures/ma_plot.png", p_ma, width = 7, height = 5, dpi = 200)

# 4d. volcano plot
vd <- tab[!is.na(tab$padj) & !is.na(tab$log2FC_shrunk), ]
vd$nlp <- -log10(pmax(vd$padj, 1e-300))
n_half <- ceiling(N_LABEL_VOLCANO / 2)
lab_genes <- rbind(head(vd[vd$direction == up_lab, ], n_half),
                   head(vd[vd$direction == down_lab, ], n_half))
p_vol <- ggplot(vd, aes(log2FC_shrunk, nlp, colour = direction)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab_genes, aes(label = label), colour = "black",
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  labs(x = "log2 fold change (shrunken)", y = "-log10 adjusted p-value", colour = NULL,
       title = paste0(GSE, ": ", GROUP_TEST, " vs ", GROUP_REF, " (paired design)")) +
  theme_minimal(base_size = 11)
ggsave("figures/volcano_plot.png", p_vol, width = 7.5, height = 6, dpi = 200)

# 4e. heatmap of the most significant genes (VST, each gene centred and scaled)
top <- head(sig_tab$gene_id, N_HEATMAP_GENES)
if (length(top) >= 2) {
  hm <- vst[top, , drop = FALSE]
  rownames(hm) <- make.unique(tab$label[match(top, tab$gene_id)])
  hm <- t(scale(t(hm)))
  ann_col <- data.frame(group = meta[colnames(hm), "group"], row.names = colnames(hm))
  ann_colors <- list(group = grp_cols)
  pheatmap(hm, annotation_col = ann_col, annotation_colors = ann_colors,
           show_colnames = FALSE, fontsize_row = 7,
           color = colorRampPalette(c("#0072B2", "white", "#D55E00"))(100),
           main = paste0("Top ", length(top), " genes by adjusted p-value"),
           filename = "figures/top_genes_heatmap.png", width = 8, height = 9)
}

# 4f. paired view of the top genes: one line per patient
top6 <- head(sig_tab$gene_id, 6)
if (length(top6) > 0) {
  long <- do.call(rbind, lapply(top6, function(g) {
    data.frame(gene = tab$label[match(g, tab$gene_id)], sample_id = colnames(nc),
               expr = log2(nc[g, ] + 1), stringsAsFactors = FALSE)
  }))
  long$tissue  <- meta$group[match(long$sample_id, meta$sample_id)]
  long$patient <- meta$patient[match(long$sample_id, meta$sample_id)]
  long$gene    <- factor(long$gene, levels = unique(long$gene))
  p_pair <- ggplot(long, aes(tissue, expr)) +
    geom_line(aes(group = patient), colour = "grey70") +
    geom_point(aes(colour = tissue), size = 2) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_colour_manual(values = grp_cols) +
    labs(x = NULL, y = "log2(normalised count + 1)", colour = NULL,
         title = "Top genes: each line is one patient") +
    theme_minimal(base_size = 11)
  ggsave("figures/top_genes_paired.png", p_pair, width = 8, height = 5.5, dpi = 200)
}

# ------------------------------------------------ 5. sensitivity analysis
excl <- SENS_EXCLUDE_PATIENT
if (excl %in% as.character(meta$patient)) {
  keep  <- as.character(meta$patient) != excl
  fit_s <- run_deseq(mat[, keep, drop = FALSE], meta[colnames(mat)[keep], ])
  tab_s <- make_table(fit_s)
  write.csv(tab_s, paste0("results/deseq2_sensitivity_without_", excl, ".csv"), row.names = FALSE)

  m <- merge(tab[, c("gene_id", "log2FC_shrunk", "direction")],
             tab_s[, c("gene_id", "log2FC_shrunk", "direction")],
             by = "gene_id", suffixes = c("_all", "_sens"))
  r <- cor(m$log2FC_shrunk_all, m$log2FC_shrunk_sens, use = "complete.obs")
  s_all  <- m$gene_id[m$direction_all  != "Not significant"]
  s_sens <- m$gene_id[m$direction_sens != "Not significant"]
  both   <- intersect(s_all, s_sens)
  jacc   <- if (length(union(s_all, s_sens)) > 0) length(both) / length(union(s_all, s_sens)) else NA

  writeLines(c(
    paste0("# Sensitivity analysis: without patient ", excl, " (both samples)"), "",
    paste0("- Samples: ", sum(keep), " (all: ", ncol(mat), ")"),
    paste0("- Correlation of shrunken log2FC between both analyses: ", sprintf("%.3f", r)),
    paste0("- Significant genes, all samples: ", length(s_all), "; without ", excl, ": ", length(s_sens)),
    paste0("- Significant in both: ", length(both), " (overlap = ",
           ifelse(is.na(jacc), "NA", sprintf("%.1f%%", 100 * jacc)), " of the union)"),
    "", "Interpretation: a high correlation and a large overlap mean that the conclusions",
    "do not depend on this patient."
  ), "results/sensitivity_summary.md")

  p_s <- ggplot(m, aes(log2FC_shrunk_all, log2FC_shrunk_sens)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, colour = "#D55E00", linetype = "dashed") +
    labs(x = "log2FC, all samples", y = paste0("log2FC, without ", excl),
         title = paste0("Sensitivity analysis (r = ", sprintf("%.3f", r), ")")) +
    theme_minimal(base_size = 11)
  ggsave("figures/sensitivity_lfc_scatter.png", p_s, width = 5.5, height = 5, dpi = 200)
} else {
  message("Patient ", excl, " not found: sensitivity analysis skipped.")
}

writeLines(capture.output(sessionInfo()), "metadata/sessionInfo_deseq2.txt")
message(paste(summ, collapse = "\n"))
