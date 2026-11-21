# Public RNA-seq Disease Response Explorer

**Author:** Malek Aloui | DrugIT bioinformatics mini-project (2026)

Code: https://github.com/MALEK9639/drugit-rnaseq

## Biological question

Which genes and pathways differ between normal colon tissue and primary colorectal cancer in a public bulk RNA-seq dataset (GEO: GSE50760, 18 patients, paired normal / tumour)?

## Status

* \[x] Week 1: setup, dataset selection (see `docs/` and `metadata/dataset\_comparison.md`)
* \[x] Week 2: scripts 01-03 (import, cleaning, quality control)
* \[x] Week 3: scripts 04-05 (DESeq2, volcano / MA / heatmap, pathway enrichment)
* \[x] Week 4: Shiny app (deployed), report and presentation (see `report/`)Installation
1. Install R (>= 4.3) and RStudio.
2. From the project root: `source("scripts/00\_install\_packages.R")` then `source("scripts/00b\_install\_week3\_packages.R")`

## Execution order

Run from the project root (the folder containing this README).

|Order|Script|Output|
|-|-|-|
|0|`scripts/00\_config.R` (sourced by the others)|settings|
|1|`scripts/01\_download\_or\_import\_data.R`|`data/raw/`, `metadata/retrieval\_log.csv`|
|2|`scripts/02\_clean\_metadata.R`|`data/processed/counts\_clean.csv`, `metadata/sample\_metadata\_clean.csv`, `metadata/cleaning\_log.txt`|
|2b|`scripts/02b\_add\_gene\_symbols.R`|`metadata/gene\_annotation.csv` (gene symbols, display only)|
|3|`scripts/03\_quality\_control.R`|`figures/library\_size.png`, `figures/sample\_distance\_heatmap.png`, `figures/pca\_plot.png`, `results/week2\_qc\_summary.md`|

| 4 | `scripts/04\_differential\_expression.R` | `results/deseq2\_all\_genes.csv`, `results/deseq2\_significant\_genes.csv`, volcano / MA / heatmap / paired plots, sensitivity analysis without P17 |
| 5 | `scripts/05\_pathway\_enrichment.R` | `results/go\_bp\_\*.csv`, `results/gsea\_go\_bp.csv`, enrichment figures |

| 6 | `scripts/06\_prepare\_app\_data.R` | `app/data/\*.rds` (compact data for the app) |
| App | `shiny::runApp("app")` | interactive explorer (PCA, volcano, gene explorer, enrichment) |

## Data

Raw count matrices are not committed (`.gitignore`); script 01 re-downloads them and logs accession, URL and retrieval date.

## Main findings

Paired analysis of 18 patients (normal colon vs primary tumour, DESeq2, design `\~ patient + group`), 18,894 genes tested.

* **2,211 genes** pass adjusted p < 0.05 and |shrunken log2FC| >= 1: **646 up** and **1,565 down** in the tumour.
* **Up in tumour:** KRT80, INHBA, ETV4, CST1, COL11A1, CLDN1, CDH3, ESM1, MMP7, COL10A1.
* **Down in tumour:** BEST4, OTOP2, CA7, SPIB, OTOP3, UNC5C, SH2D6, ABCG2, CLCA4-AS1, MYOC (markers of mature colonic epithelium).
* **Over-representation (GO BP):** up genes are enriched for extracellular matrix organisation, angiogenesis and epithelial-mesenchymal transition; down genes for immunoglobulin-mediated immunity, ion homeostasis and sodium transport.
* **GSEA (all genes ranked):** ribosome biogenesis, rRNA processing, DNA replication and cell-cycle processes are higher in the tumour (proliferation); immune, ion-transport and lipid-metabolism processes are lower. Proliferation genes change moderately, so GSEA detects them while thresholded ORA does not.
* **Sensitivity analysis:** removing patient P17 (flagged as an outlier in QC) gives a log2FC correlation of 0.997; 2,159 of the 2,211 genes remain significant (87% overlap of the union).

## Interactive app

Run `shiny::runApp("app")` after `source("scripts/06\_prepare\_app\_data.R")`. Tabs: overview, PCA, volcano with adjustable thresholds, per-gene paired plot, enrichment.
Live app: https://malekaloui.shinyapps.io/rnaseq-disease-response-explorer/

## Limitations

* Bulk tissue: part of the signal reflects differences in cell composition (for example fewer plasma cells and immune cells in the tumour samples), not only changes inside tumour cells.
* All patients have stage IV disease with liver metastasis, so results may not generalise to earlier stages.
* 18 patients: moderate power for small effects; thresholds (padj 0.05, |log2FC| 1) are conventions.
* Associations only; no experimental validation. Pathway results depend on GO annotations.
* Three samples were flagged in QC (GSM1228203, GSM1228212, GSM1228213); none was removed, and the P17 sensitivity analysis shows the conclusions are stable.

## Reference

Kim SK et al. (2014). A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. GEO: GSE50760.

