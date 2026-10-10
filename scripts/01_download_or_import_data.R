# 01_download_or_import_data.R
# Goal: download the raw count matrix + sample metadata for the approved dataset,
#       keep untouched copies in data/raw, and log accession / URL / retrieval date.
# Input : none (internet access)
# Output: data/raw/<GSE>_raw_counts_<genome>_NCBI.tsv.gz
#         data/raw/Human.<genome>.annot.tsv.gz
#         data/raw/geo_sample_metadata_raw.csv
#         metadata/retrieval_log.csv

source("scripts/00_config.R")
suppressPackageStartupMessages(library(GEOquery))
options(timeout = 900)

counts_name <- sprintf("%s_raw_counts_%s_NCBI.tsv.gz", GSE, GENOME_TAG)
annot_name  <- sprintf("Human.%s.annot.tsv.gz", GENOME_TAG)

# NCBI-generated raw counts (same pipeline for every human GEO RNA-seq series).
counts_url <- sprintf(
  "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s",
  GSE, counts_name)
annot_url <- sprintf(
  "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=%s",
  annot_name)

fetch <- function(url, dest) {
  if (file.exists(dest)) { message("Already present, not re-downloaded: ", dest); return(invisible(TRUE)) }
  ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = FALSE); TRUE },
                 error = function(e) { message("Download failed: ", conditionMessage(e)); FALSE },
                 warning = function(w) { message("Download warning: ", conditionMessage(w)); FALSE })
  if (!ok && file.exists(dest)) file.remove(dest)
  if (!ok) message("--> Download the file manually from the GEO page of ", GSE,
                   " ('Download' section, 'RNA-seq raw counts matrix') and put it in data/raw/")
  invisible(ok)
}

fetch(counts_url, file.path("data/raw", counts_name))
fetch(annot_url,  file.path("data/raw", annot_name))

# Sample metadata from the GEO series matrix
message("Fetching sample metadata for ", GSE, " ...")
gse <- getGEO(GSE, GSEMatrix = TRUE, getGPL = FALSE)
if (is.list(gse)) gse <- gse[[1]]
pd <- Biobase::pData(gse)
write.csv(pd, "data/raw/geo_sample_metadata_raw.csv")
message("Samples in metadata: ", nrow(pd))

# Retrieval log (reproducibility requirement)
log_df <- data.frame(
  accession      = GSE,
  file           = c(counts_name, annot_name, "geo_sample_metadata_raw.csv"),
  source_url     = c(counts_url, annot_url, paste0("GEOquery::getGEO('", GSE, "')")),
  retrieval_date = as.character(Sys.Date()),
  R_version      = R.version.string,
  GEOquery       = as.character(packageVersion("GEOquery")),
  stringsAsFactors = FALSE
)
write.csv(log_df, "metadata/retrieval_log.csv", row.names = FALSE)
message("Done. See metadata/retrieval_log.csv")
