# 02_clean_metadata.R
# Goal: build a clean sample table and a clean count matrix (data/processed),
#       and write every check to metadata/cleaning_log.txt.
# Input : data/raw/* from script 01
# Output: metadata/sample_metadata_clean.csv
#         data/processed/counts_clean.csv
#         metadata/gene_annotation.csv
#         metadata/cleaning_log.txt

source("scripts/00_config.R")

logf <- "metadata/cleaning_log.txt"
log_lines <- character()
say <- function(...) { msg <- paste0(...); message(msg); log_lines <<- c(log_lines, msg) }

# ------------------------------------------------------------------ 1. metadata
pd <- read.csv("data/raw/geo_sample_metadata_raw.csv", row.names = 1,
               check.names = FALSE, stringsAsFactors = FALSE)
say("Samples in GEO metadata: ", nrow(pd))
say("Metadata columns: ", paste(colnames(pd), collapse = " | "))

# >>> INSPECT FIRST: look at these before trusting the automatic labelling below <<<
print(head(pd[, intersect(c("title", "source_name_ch1"), colnames(pd)), drop = FALSE], 10))
chr_cols <- grep("^characteristics", colnames(pd))
if (length(chr_cols)) print(head(pd[, chr_cols, drop = FALSE], 5))

# Tissue class, derived from title + source name + characteristics (lower case).
# If the printed table below is wrong or contains NA, EDIT these rules.
txt <- tolower(apply(pd[, intersect(c("title", "source_name_ch1", colnames(pd)[chr_cols]),
                                    colnames(pd)), drop = FALSE], 1, paste, collapse = " "))
pd$tissue_class <- ifelse(grepl("normal", txt), "normal",
                   ifelse(grepl("metasta", txt), "metastasis",
                   ifelse(grepl("primary|tumou?r|cancer|crc", txt), "tumor", NA)))
say("Tissue classes found:")
print(table(pd$tissue_class, useNA = "ifany"))
if (anyNA(pd$tissue_class))
  stop("Some samples have no tissue class: edit the rules in section 1 after inspecting the metadata.")

# Patient identifier (needed for the paired design). Taken from the first number in the title.
# VERIFY on the printout: each patient should have one sample per tissue class.
pd$patient <- paste0("P", sub(".*?([0-9]+).*", "\\1", pd$title, perl = TRUE))
say("Patients detected: ", length(unique(pd$patient)))
print(table(pd$patient, pd$tissue_class))

# Keep only the two groups of the main comparison (metastasis is out of scope).
meta <- pd[pd$tissue_class %in% c(GROUP_REF, GROUP_TEST), , drop = FALSE]
meta$sample_id <- rownames(meta)                       # GSM accession
meta$group     <- factor(meta$tissue_class, levels = c(GROUP_REF, GROUP_TEST))
meta <- meta[, c("sample_id", "title", "patient", "group")]

# Checks required by the assignment
say("Samples kept for the main comparison: ", nrow(meta))
say("Group sizes: ", paste(names(table(meta$group)), table(meta$group), collapse = ", "))
say("Missing group labels: ", sum(is.na(meta$group)))
say("Duplicated sample IDs: ", sum(duplicated(meta$sample_id)))
say("Patients missing one of the two groups: ",
    sum(rowSums(table(meta$patient, meta$group) > 0) < 2))
stopifnot(!anyNA(meta$group), !any(duplicated(meta$sample_id)))
write.csv(meta, "metadata/sample_metadata_clean.csv", row.names = FALSE)

# -------------------------------------------------------------------- 2. counts
counts_file <- sprintf("data/raw/%s_raw_counts_%s_NCBI.tsv.gz", GSE, GENOME_TAG)
cnt <- read.delim(gzfile(counts_file), check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- as.character(cnt[[1]])
mat <- as.matrix(cnt[, -1, drop = FALSE])
say("Raw matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

# Non-gene summary rows (e.g. HTSeq counters __no_feature). NCBI matrices use numeric GeneIDs.
non_gene <- grepl("^__|^N_|no_feature|ambiguous", gene_ids)
say("Non-gene summary rows removed: ", sum(non_gene))
mat <- mat[!non_gene, , drop = FALSE]; gene_ids <- gene_ids[!non_gene]

say("Duplicated gene identifiers: ", sum(duplicated(gene_ids)))
stopifnot(!any(duplicated(gene_ids)))
rownames(mat) <- gene_ids

say("Missing values in counts: ", sum(is.na(mat)))
stopifnot(!anyNA(mat))
say("Negative values: ", sum(mat < 0))
say("Non-integer values: ", sum(mat != round(mat)))
mat <- round(mat)                                      # rounding only if needed; logged above
storage.mode(mat) <- "integer"

# Sample matching between count matrix and metadata
say("Samples in count matrix but not in metadata (expected: the metastasis samples): ",
    length(setdiff(colnames(mat), meta$sample_id)))
missing_in_counts <- setdiff(meta$sample_id, colnames(mat))
say("Metadata samples missing from the count matrix: ", length(missing_in_counts))
stopifnot(length(missing_in_counts) == 0)
mat <- mat[, meta$sample_id, drop = FALSE]
stopifnot(identical(colnames(mat), meta$sample_id))
say("Sample names match metadata: TRUE")
write.csv(mat, "data/processed/counts_clean.csv")

# ----------------------------------------------------------- 3. gene identifiers
# NCBI matrices already use Entrez GeneIDs; we do NOT convert them for the analysis.
# Symbols are added only for display, from the NCBI annotation table.
annot_file <- sprintf("data/raw/Human.%s.annot.tsv.gz", GENOME_TAG)
if (file.exists(annot_file)) {
  ann <- read.delim(gzfile(annot_file), stringsAsFactors = FALSE, quote = "")
  id_col  <- if ("GeneID" %in% names(ann)) "GeneID" else names(ann)[1]
  sym_col <- if ("Symbol" %in% names(ann)) "Symbol" else NA
  ann_out <- data.frame(gene_id = as.character(ann[[id_col]]),
                        symbol  = if (!is.na(sym_col)) ann[[sym_col]] else NA,
                        stringsAsFactors = FALSE)
  ann_out <- ann_out[!duplicated(ann_out$gene_id), ]
  say("Genes without a symbol in the annotation: ",
      sum(!rownames(mat) %in% ann_out$gene_id | is.na(ann_out$symbol[match(rownames(mat), ann_out$gene_id)])))
  write.csv(ann_out, "metadata/gene_annotation.csv", row.names = FALSE)
  say("ID method: Entrez GeneID kept; symbols added from ", basename(annot_file), " (display only).")
} else say("Annotation file not found: symbols will be missing.")

writeLines(log_lines, logf)
message("Done. See ", logf)
