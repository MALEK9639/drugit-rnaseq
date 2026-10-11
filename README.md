# Public RNA-seq Disease Response Explorer

**Author:** Malek Aloui | DrugIT bioinformatics mini-project (2026)

Code: https://github.com/MALEK9639/drugit-rnaseq

## Biological question
Which genes and pathways differ between normal colon tissue and primary colorectal cancer in a public bulk RNA-seq dataset (GEO: GSE50760, 18 patients, paired normal / tumour)?

## Status
- [x] Week 1: setup, dataset selection (see `docs/` and `metadata/dataset_comparison.md`)
- [x] Week 2: scripts 01-03 (import, cleaning, quality control)
- [x] Week 3: scripts 04-05 (DESeq2, volcano / MA / heatmap, pathway enrichment)
- [ ] Week 4: Shiny app (in progress), report, presentation

## Installation
1. Install R (>= 4.3) and RStudio.
2. From the project root: `source("scripts/00_install_packages.R")` then `source("scripts/00b_install_week3_packages.R")`

## Execution order
Run from the project root (the folder containing this README).

| Order | Script | Output |
|---|---|---|
| 0 | `scripts/00_config.R` (sourced by the others) | settings |
| 1 | `scripts/01_download_or_import_data.R` | `data/raw/`, `metadata/retrieval_log.csv` |
| 2 | `scripts/02_clean_metadata.R` | `data/processed/counts_clean.csv`, `metadata/sample_metadata_clean.csv`, `metadata/cleaning_log.txt` |
| 2b | `scripts/02b_add_gene_symbols.R` | `metadata/gene_annotation.csv` (gene symbols, display only) |
| 3 | `scripts/03_quality_control.R` | `figures/library_size.png`, `figures/sample_distance_heatmap.png`, `figures/pca_plot.png`, `results/week2_qc_summary.md` |

| 4 | `scripts/04_differential_expression.R` | `results/deseq2_all_genes.csv`, `results/deseq2_significant_genes.csv`, volcano / MA / heatmap / paired plots, sensitivity analysis without P17 |
| 5 | `scripts/05_pathway_enrichment.R` | `results/go_bp_*.csv`, `results/gsea_go_bp.csv`, enrichment figures |

| 6 | `scripts/06_prepare_app_data.R` | `app/data/*.rds` (compact data for the app) |
| App | `shiny::runApp("app")` | interactive explorer (PCA, volcano, gene explorer, enrichment) |

## Data
Raw count matrices are not committed (`.gitignore`); script 01 re-downloads them and logs accession, URL and retrieval date.

## Main findings
Paired analysis of 18 patients (normal colon vs primary tumour, DESeq2, design `~ patient + group`), 18,894 genes tested.

- **2,211 genes** pass adjusted p < 0.05 and |shrunken log2FC| >= 1: **646 up** and **1,565 down** in the tumour.
- **Up in tumour:** KRT80, INHBA, ETV4, CST1, COL11A1, CLDN1, CDH3, ESM1, MMP7, COL10A1.
- **Down in tumour:** BEST4, OTOP2, CA7, SPIB, OTOP3, UNC5C, SH2D6, ABCG2, CLCA4-AS1, MYOC (markers of mature colonic epithelium).
- **Over-representation (GO BP):** up genes are enriched for extracellular matrix organisation, angiogenesis and epithelial-mesenchymal transition; down genes for immunoglobulin-mediated immunity, ion homeostasis and sodium transport.
- **GSEA (all genes ranked):** ribosome biogenesis, rRNA processing, DNA replication and cell-cycle processes are higher in the tumour (proliferation); immune, ion-transport and lipid-metabolism processes are lower. Proliferation genes change moderately, so GSEA detects them while thresholded ORA does not.
- **Sensitivity analysis:** removing patient P17 (flagged as an outlier in QC) gives a log2FC correlation of 0.997; 2,159 of the 2,211 genes remain significant (87% overlap of the union).

## Interactive app
Run `shiny::runApp("app")` after `source("scripts/06_prepare_app_data.R")`. Tabs: overview, PCA, volcano with adjustable thresholds, per-gene paired plot, enrichment.
Live app: https://malekaloui.shinyapps.io/rnaseq-disease-response-explorer/

## Limitations
- Bulk tissue: part of the signal reflects differences in cell composition (for example fewer plasma cells and immune cells in the tumour samples), not only changes inside tumour cells.
- All patients have stage IV disease with liver metastasis, so results may not generalise to earlier stages.
- 18 patients: moderate power for small effects; thresholds (padj 0.05, |log2FC| 1) are conventions.
- Associations only; no experimental validation. Pathway results depend on GO annotations.
- Three samples were flagged in QC (GSM1228203, GSM1228212, GSM1228213); none was removed, and the P17 sensitivity analysis shows the conclusions are stable.

## Reference
Kim SK et al. (2014). A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. GEO: GSE50760.
