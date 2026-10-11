# Public RNA-seq Disease Response Explorer

**Author:** Malek Aloui | DrugIT bioinformatics mini-project (2026)

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
To be completed after week 3.

## Limitations
To be completed after week 3.
