# Dataset Selection Form

**Intern:** Malek Aloui | **Supervisor:** Koussai Salem | **Date:** ____ / ____ / 2026

| Field | Answer |
|---|---|
| Dataset accession | GSE50760 (GEO). Retrieval date and download links are recorded in `metadata/retrieval_log.csv` by script 01. |
| Biological question | Which genes and pathways differ between normal colon tissue and primary colorectal cancer tissue in the same patients? |
| Groups | `normal` (normal colon mucosa) vs `tumor` (primary colorectal cancer). Expected 18 vs 18, from 18 patients. The 18 liver-metastasis samples are excluded from the main analysis. *Confirm the exact numbers after running script 02.* |
| Count data | NCBI-generated raw count matrix `GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz` from the GEO Download page, saved untouched in `data/raw/`. *Confirm the file exists for this series.* |
| Metadata | GEO sample table (series matrix), saved as `data/raw/geo_sample_metadata_raw.csv`; cleaned version in `metadata/sample_metadata_clean.csv`. |
| Why this dataset | 1. 18 biological replicates per group, well above the minimum of 5. 2. Two clearly defined tissue groups with a patient identifier, so each tumour can be compared with its own normal tissue. 3. Cancer vs normal colon is well studied, which allows cautious interpretation of the top genes with reliable sources. |
| Main risk | All patients had stage IV disease with synchronous liver metastasis, so "primary tumour" is not representative of all colorectal cancers. The paired design also needs `~ patient + group`; I would like to confirm this is acceptable. |
| Original publication | Kim SK, Kim SY, Kim JH, Roh SA, Cho DH, et al. (2014). A nineteen gene-based risk score classifier predicts prognosis of colorectal cancer patients. *Molecular Oncology*. DOI: *to complete after checking the GEO page*. |

**Decision requested:** approval of GSE50760 (normal vs primary tumour) and of the paired design formula `~ patient + group`.
