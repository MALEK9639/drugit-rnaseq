# 03_quality_control.R
# Goal: library sizes, filtering, transformed matrix, sample-distance heatmap, PCA,
#       and an outlier FLAG table (flags only: nothing is removed without the supervisor).
# Input : data/processed/counts_clean.csv, metadata/sample_metadata_clean.csv
# Output: figures/library_size.png, figures/sample_distance_heatmap.png, figures/pca_plot.png
#         data/processed/counts_filtered.csv, data/processed/vst_matrix.csv
#         results/qc_outlier_flags.csv, results/week2_qc_summary.md

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2); library(pheatmap) })

meta <- read.csv("metadata/sample_metadata_clean.csv", stringsAsFactors = FALSE)
meta$group <- factor(meta$group, levels = c(GROUP_REF, GROUP_TEST))
mat  <- as.matrix(read.csv("data/processed/counts_clean.csv", row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(mat), meta$sample_id))
rownames(meta) <- meta$sample_id

# ---------------------------------------------------------------- library size
lib <- colSums(mat)
meta$library_size <- lib[rownames(meta)]
p_lib <- ggplot(meta, aes(x = reorder(sample_id, library_size), y = library_size / 1e6, fill = group)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Library size (million counts)", fill = "Group",
       title = paste0(GSE, ": library size per sample")) +
  theme_minimal(base_size = 9)
ggsave("figures/library_size.png", p_lib, width = 7, height = 8, dpi = 200)

# -------------------------------------------------------------------- filtering
n_min  <- min(table(meta$group))
keep   <- rowSums(mat >= MIN_COUNT) >= n_min
genes_before <- nrow(mat); genes_after <- sum(keep)
mat_f  <- mat[keep, ]
write.csv(mat_f, "data/processed/counts_filtered.csv")

# ----------------------------------------------- transformed matrix (visual only)
dds <- DESeqDataSetFromMatrix(countData = mat_f, colData = meta, design = ~ group)
vsd <- vst(dds, blind = TRUE)                       # blind: QC must not depend on the design
vmat <- assay(vsd)
write.csv(vmat, "data/processed/vst_matrix.csv")

# ----------------------------------------------------- sample-distance heatmap
d <- dist(t(vmat))
dm <- as.matrix(d)
ann_col <- data.frame(group = meta$group, row.names = rownames(meta))
pheatmap(dm, clustering_distance_rows = d, clustering_distance_cols = d,
         annotation_col = ann_col, show_rownames = FALSE, show_colnames = FALSE,
         main = paste0(GSE, ": sample-to-sample distances (VST)"),
         filename = "figures/sample_distance_heatmap.png", width = 8, height = 7)

# -------------------------------------------------------------------------- PCA
top <- order(apply(vmat, 1, var), decreasing = TRUE)[seq_len(min(N_TOP_VAR, nrow(vmat)))]
pca <- prcomp(t(vmat[top, ]), center = TRUE, scale. = FALSE)
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc  <- data.frame(pca$x[, 1:2], meta[rownames(pca$x), c("sample_id", "group", "patient")])
p_pca <- ggplot(pc, aes(PC1, PC2, colour = group, shape = group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       title = paste0(GSE, ": PCA, top ", length(top), " variable genes (VST)")) +
  theme_minimal(base_size = 11)
if (requireNamespace("ggrepel", quietly = TRUE))
  p_pca <- p_pca + ggrepel::geom_text_repel(aes(label = patient), size = 2.5, max.overlaps = 30)
ggsave("figures/pca_plot.png", p_pca, width = 7, height = 5.5, dpi = 200)

# ---------------------------------------------------------------- outlier flags
mad_z <- function(x) abs(x - median(x)) / (mad(x) + 1e-9)
cent  <- aggregate(cbind(PC1, PC2) ~ group, pc, mean)
pc$dist_centroid <- mapply(function(g, x, y) {
  c0 <- cent[as.character(cent$group) == g, ]; sqrt((x - c0$PC1)^2 + (y - c0$PC2)^2)
}, as.character(pc$group), pc$PC1, pc$PC2)
# Distance to the centroid of the sample's OWN group, scored within each group
# (a pooled score would hide outliers because tumours are naturally more spread out).
d_c   <- pc$dist_centroid[match(rownames(meta), rownames(pc))]
pca_z <- ave(d_c, as.character(meta$group), FUN = mad_z)
flags <- data.frame(
  sample_id       = rownames(meta),
  group           = meta$group,
  patient         = meta$patient,
  library_size    = meta$library_size,
  lib_below_half_median = meta$library_size < 0.5 * median(meta$library_size),
  pca_dist_z      = round(pca_z, 1),
  pca_dist_flag   = pca_z > 3.5
)
flags$flagged <- flags$lib_below_half_median | flags$pca_dist_flag
write.csv(flags, "results/qc_outlier_flags.csv", row.names = FALSE)

# ------------------------------------------------------------- summary for email
summ <- c(
  paste0("# Week 2 QC summary (", GSE, ")"), "",
  paste0("- Samples: ", nrow(meta), " (", paste(names(table(meta$group)), table(meta$group), collapse = ", "), ")"),
  paste0("- Missing values in count matrix: ", sum(is.na(mat))),
  paste0("- Filtering rule: count >= ", MIN_COUNT, " in at least ", n_min, " samples (size of the smallest group)"),
  paste0("- Genes before / after filtering: ", genes_before, " / ", genes_after),
  paste0("- Library size: min ", round(min(lib) / 1e6, 1), "M, median ", round(median(lib) / 1e6, 1),
         "M, max ", round(max(lib) / 1e6, 1), "M"),
  paste0("- PCA: PC1 ", pve[1], "%, PC2 ", pve[2], "%"),
  paste0("- Samples flagged as possible outliers: ", sum(flags$flagged),
         if (any(flags$flagged)) paste0(" (", paste(flags$sample_id[flags$flagged], collapse = ", "), ")") else ""),
  "- No sample is removed. Any exclusion must be justified in writing and discussed with the supervisor."
)
writeLines(summ, "results/week2_qc_summary.md")
writeLines(capture.output(sessionInfo()), "metadata/sessionInfo_qc.txt")
message(paste(summ, collapse = "\n"))
