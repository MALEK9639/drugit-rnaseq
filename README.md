# Public RNA-seq Disease Response Explorer

**Intern:** Malek Aloui | **Organisation:** DrugIT | **Supervisor:** Koussai Salem
**Internship dates:** 15 / 9 / 2026 to 15 / 10 / 2026

## Biological question

Which genes and pathways differ between normal colon tissue and primary colorectal cancer in a public bulk RNA-seq dataset (GEO: GSE50760, pending supervisor approval)?

## Status

* \[x] Week 1: setup, dataset selection (see `docs/` and `metadata/dataset\_comparison.md`)
* \[x] Week 2: scripts 01-03 (import, cleaning, quality control)
* \[ ] Week 3: DESeq2, volcano / MA / heatmap, pathway enrichment
* \[ ] Week 4: Shiny app, report, presentation

## Installation

1. Install R (>= 4.3) and RStudio.
2. From the project root: `source("scripts/00\_install\_packages.R")`

## Execution order

Run from the project root (the folder containing this README).

|Order|Script|Output|
|-|-|-|
|0|`scripts/00\_config.R` (sourced by the others)|settings|
|1|`scripts/01\_download\_or\_import\_data.R`|`data/raw/`, `metadata/retrieval\_log.csv`|
|2|`scripts/02\_clean\_metadata.R`|`data/processed/counts\_clean.csv`, `metadata/sample\_metadata\_clean.csv`, `metadata/cleaning\_log.txt`|
|3|`scripts/03\_quality\_control.R`|`figures/library\_size.png`, `figures/sample\_distance\_heatmap.png`, `figures/pca\_plot.png`, `results/week2\_qc\_summary.md`|

Scripts 04-06 and the Shiny app will be added in weeks 3 and 4.

## Data

Raw count matrices are not committed (`.gitignore`); script 01 re-downloads them and logs accession, URL and retrieval date.

## Main findings

To be completed after week 3.

## Limitations

To be completed after week 3.

