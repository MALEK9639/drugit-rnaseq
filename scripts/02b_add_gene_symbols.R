# 02b_add_gene_symbols.R
# Goal: add gene symbols and names (display only) to the Entrez GeneIDs of the count matrix,
#       using the Bioconductor package org.Hs.eg.db (no download needed).
# Input : data/processed/counts_clean.csv (from script 02)
# Output: metadata/gene_annotation.csv (columns: gene_id, symbol, gene_name)

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

ids <- rownames(read.csv("data/processed/counts_clean.csv", row.names = 1, check.names = FALSE))

symbol    <- mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "SYMBOL",   multiVals = "first")
gene_name <- mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "GENENAME", multiVals = "first")

ann <- data.frame(gene_id = ids, symbol = unname(symbol), gene_name = unname(gene_name),
                  stringsAsFactors = FALSE)
write.csv(ann, "metadata/gene_annotation.csv", row.names = FALSE)

message("Genes: ", nrow(ann))
message("With a symbol: ", sum(!is.na(ann$symbol)))
message("Without a symbol (will be shown by their GeneID): ", sum(is.na(ann$symbol)))