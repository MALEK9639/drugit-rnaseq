# 00b_install_week3_packages.R -- run once, before script 04.
# apeglm: better (shrunken) log2 fold-change estimates for low-count genes.
# bslib: nicer layout for the Shiny app (week 4).
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("apeglm", ask = FALSE, update = FALSE)
install.packages("bslib")

# Check that everything needed for weeks 3-4 can be loaded
pkgs <- c("DESeq2", "apeglm", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
          "ggplot2", "ggrepel", "pheatmap", "shiny", "DT", "plotly", "bslib")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
print(ok)
if (!all(ok)) message("Missing packages: ", paste(names(ok)[!ok], collapse = ", "))
