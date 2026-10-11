# 06_prepare_app_data.R
# Goal: build compact data files for the Shiny app, so that the app/ folder is self-contained
#       (it can be deployed without the raw data or the other scripts).
# Input : results/deseq2_all_genes.csv, results/go_bp_*.csv, results/gsea_go_bp.csv,
#         data/processed/vst_matrix.csv, metadata/sample_metadata_clean.csv
# Output: app/data/results.rds, metadata.rds, expression_vst.rds, pca.rds, enrichment.rds, params.rds

source("scripts/00_config.R")
out <- "app/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

up_lab   <- paste("Up in", GROUP_TEST)
down_lab <- paste("Down in", GROUP_TEST)

# ------------------------------------------------ 1. differential expression table
res <- read.csv("results/deseq2_all_genes.csv", colClasses = c(gene_id = "character"))
need <- c("gene_id", "label", "symbol", "gene_name", "baseMean", "log2FC",
          "log2FC_shrunk", "pvalue", "padj", "stat", "direction")
stopifnot(all(need %in% names(res)))
res <- res[, need]
res$gene_name[is.na(res$gene_name)] <- ""
res$label[is.na(res$label) | res$label == ""] <- res$gene_id[is.na(res$label) | res$label == ""]

# ------------------------------------------------ 2. metadata and expression
meta <- read.csv("metadata/sample_metadata_clean.csv", stringsAsFactors = FALSE)
meta <- meta[, c("sample_id", "patient", "group")]
meta$group <- factor(meta$group, levels = c(GROUP_REF, GROUP_TEST))

vst <- as.matrix(read.csv("data/processed/vst_matrix.csv", row.names = 1, check.names = FALSE))
stopifnot(all(meta$sample_id %in% colnames(vst)))
vst <- vst[, meta$sample_id, drop = FALSE]
vst <- vst[rownames(vst) %in% res$gene_id, , drop = FALSE]
vst <- round(vst, 3)

# ------------------------------------------------ 3. PCA (top variable genes)
v   <- apply(vst, 1, var)
top <- names(sort(v, decreasing = TRUE))[seq_len(min(N_TOP_VAR, length(v)))]
pc  <- prcomp(t(vst[top, ]), center = TRUE, scale. = FALSE)
npc <- min(5, ncol(pc$x))
scores <- data.frame(sample_id = rownames(pc$x), pc$x[, seq_len(npc), drop = FALSE],
                     stringsAsFactors = FALSE)
rownames(scores) <- NULL
scores$patient <- meta$patient[match(scores$sample_id, meta$sample_id)]
scores$group   <- meta$group[match(scores$sample_id, meta$sample_id)]
var_expl <- 100 * pc$sdev[seq_len(npc)]^2 / sum(pc$sdev^2)
names(var_expl) <- paste0("PC", seq_len(npc))
pca <- list(scores = scores, var = var_expl)

# ------------------------------------------------ 4. enrichment tables
read_enr <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
add_ratio <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(NULL)
  parts <- strsplit(d$GeneRatio, "/")
  d$GeneRatioNum <- vapply(parts, function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  d
}
gsea <- read_enr("results/gsea_go_bp.csv")
if (!is.null(gsea) && nrow(gsea) == 0) gsea <- NULL
enr <- list(up   = add_ratio(read_enr("results/go_bp_up_in_tumor.csv")),
            down = add_ratio(read_enr("results/go_bp_down_in_tumor.csv")),
            gsea = gsea)

# ------------------------------------------------ 5. parameters shown in the app
prm <- list(
  gse = GSE, group_ref = GROUP_REF, group_test = GROUP_TEST,
  padj_cutoff = PADJ_CUTOFF, lfc_cutoff = LFC_CUTOFF,
  n_samples = nrow(meta), n_patients = length(unique(meta$patient)),
  n_genes = sum(!is.na(res$padj)),
  n_up = sum(res$direction == up_lab, na.rm = TRUE),
  n_down = sum(res$direction == down_lab, na.rm = TRUE)
)

saveRDS(res,  file.path(out, "results.rds"),        compress = "xz")
saveRDS(meta, file.path(out, "metadata.rds"),       compress = "xz")
saveRDS(vst,  file.path(out, "expression_vst.rds"), compress = "xz")
saveRDS(pca,  file.path(out, "pca.rds"),            compress = "xz")
saveRDS(enr,  file.path(out, "enrichment.rds"),     compress = "xz")
saveRDS(prm,  file.path(out, "params.rds"),         compress = "xz")

message("App data written to ", out, ":")
print(file.info(list.files(out, full.names = TRUE))[, "size", drop = FALSE])
