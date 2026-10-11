# 05_pathway_enrichment.R
# Goal: which biological processes are over-represented among the differentially expressed
#       genes (over-representation analysis, ORA) and which are shifted across the whole
#       ranked gene list (GSEA). Gene Ontology, Biological Process.
# Input : results/deseq2_all_genes.csv (from script 04)
# Output: results/go_bp_<tag>.csv, results/gsea_go_bp.csv, results/enrichment_summary.md,
#         figures/go_bp_<tag>_dotplot.png, figures/gsea_go_dotplot.png, figures/gsea_top3.png

source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)
})
set.seed(1)

tab <- read.csv("results/deseq2_all_genes.csv", colClasses = c(gene_id = "character"))
up_lab   <- paste("Up in", GROUP_TEST)
down_lab <- paste("Down in", GROUP_TEST)

universe <- tab$gene_id[!is.na(tab$padj)]          # every gene that could have been detected
up   <- tab$gene_id[tab$direction == up_lab]
down <- tab$gene_id[tab$direction == down_lab]
message("Universe: ", length(universe), " | up: ", length(up), " | down: ", length(down))

# ------------------------------------------------ 1. over-representation (ORA)
run_ora <- function(genes, tag, title) {
  if (length(genes) < 10) { message("Too few genes for ", tag, ": skipped."); return(NULL) }
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No enriched GO term for ", tag); return(NULL)
  }
  # merge near-identical GO terms so the figure is readable
  ego <- tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
                  error = function(e) ego)
  write.csv(as.data.frame(ego), sprintf("results/go_bp_%s.csv", tag), row.names = FALSE)
  p <- dotplot(ego, showCategory = 15) + ggtitle(title)
  ggsave(sprintf("figures/go_bp_%s_dotplot.png", tag), p, width = 9, height = 7, dpi = 200)
  ego
}
ego_up   <- run_ora(up,   "up_in_tumor",   paste0("GO Biological Process: ", up_lab))
ego_down <- run_ora(down, "down_in_tumor", paste0("GO Biological Process: ", down_lab))

# ----------------------------------------------------------- 2. GSEA
# Rank ALL tested genes by the Wald statistic (sign = direction, size = confidence).
gl <- tab$stat
names(gl) <- tab$gene_id
gl <- sort(gl[!is.na(gl)], decreasing = TRUE)

gse <- tryCatch(
  gseGO(geneList = gl, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
        minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
        eps = 0, seed = TRUE, verbose = FALSE),
  error = function(e) { message("GSEA failed: ", conditionMessage(e)); NULL })

if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
  gse <- tryCatch(simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                  error = function(e) gse)
  gse <- tryCatch(setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                  error = function(e) gse)
  write.csv(as.data.frame(gse), "results/gsea_go_bp.csv", row.names = FALSE)

  tryCatch({
    p_g <- dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
      ggtitle(paste0("GSEA, GO Biological Process (positive NES = higher in ", GROUP_TEST, ")"))
    ggsave("figures/gsea_go_dotplot.png", p_g, width = 11, height = 8, dpi = 200)
  }, error = function(e) message("GSEA dotplot failed: ", conditionMessage(e)))

  tryCatch({
    p_top <- gseaplot2(gse, geneSetID = seq_len(min(3, nrow(as.data.frame(gse)))))
    ggsave("figures/gsea_top3.png", p_top, width = 9, height = 6, dpi = 200)
  }, error = function(e) message("gseaplot2 failed: ", conditionMessage(e)))
} else {
  message("No significant GSEA term.")
}

# ----------------------------------------------------------- 3. summary
top_terms <- function(x, n = 8) {
  if (is.null(x)) return("- none")
  d <- as.data.frame(x)
  d <- d[order(d$p.adjust), ]
  d <- head(d, n)
  paste0("- ", d$Description, " (padj ", sprintf("%.1e", d$p.adjust), ")")
}
gse_pos <- gse_neg <- NULL
if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
  gd <- as.data.frame(gse)
  gse_pos <- gd[gd$NES > 0, ]; gse_neg <- gd[gd$NES < 0, ]
}
writeLines(c(
  "# Pathway enrichment summary (GO Biological Process)", "",
  paste0("- Universe: ", length(universe), " genes; ", up_lab, ": ", length(up),
         "; ", down_lab, ": ", length(down)), "",
  paste0("## ORA, ", tolower(up_lab)), top_terms(ego_up), "",
  paste0("## ORA, ", tolower(down_lab)), top_terms(ego_down), "",
  paste0("## GSEA, higher in ", GROUP_TEST, " (NES > 0)"),
  if (is.null(gse_pos) || nrow(gse_pos) == 0) "- none" else top_terms(gse_pos), "",
  paste0("## GSEA, lower in ", GROUP_TEST, " (NES < 0)"),
  if (is.null(gse_neg) || nrow(gse_neg) == 0) "- none" else top_terms(gse_neg)
), "results/enrichment_summary.md")

writeLines(capture.output(sessionInfo()), "metadata/sessionInfo_enrichment.txt")
message(paste(readLines("results/enrichment_summary.md"), collapse = "\n"))
