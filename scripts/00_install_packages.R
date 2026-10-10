# 00_install_packages.R -- run once. Installs every package used in the project.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(
  c("DESeq2", "GEOquery", "clusterProfiler", "org.Hs.eg.db", "enrichplot"),
  ask = FALSE, update = FALSE
)
install.packages(c("ggplot2", "ggrepel", "pheatmap", "shiny", "DT", "plotly"))

# Record the versions used (paste into the report, Methods section)
writeLines(capture.output(sessionInfo()), "metadata/sessionInfo_install.txt")
