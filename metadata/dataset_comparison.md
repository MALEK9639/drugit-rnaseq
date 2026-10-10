# Three-dataset comparison

Status of each cell: **[checked]** = seen in a source while preparing this table; **[TO VERIFY]** = open the GEO page and confirm before sending (GEO pages were not reachable by the automatic tool, so some fields come from papers or mirrors).

| Criterion | A. GSE50760 | B. GSE147507 (NHBE subset) | C. GSE183947 |
|---|---|---|---|
| Organism | Homo sapiens [checked] | Homo sapiens [checked] | Homo sapiens (per citing paper) |
| Tissue / cells | Colon tissue, colorectal cancer patients | Primary human bronchial epithelial cells (NHBE) | Breast tissue |
| Comparison | Normal colon vs primary colorectal cancer | Mock-treated vs SARS-CoV-2 infected | Normal vs malignant breast tissue |
| Design type | Disease vs healthy control (paired by patient) | Infected vs uninfected | Disease vs healthy control (paired) |
| Samples | 54 in the series (18 patients x normal / primary / liver metastasis) [checked]; 18 normal + 18 primary used | 3 mock vs 3 infected for NHBE (Series 1) [checked] | About 60 samples, i.e. matched normal and malignant pairs (per citing paper) [TO VERIFY] |
| Meets 5 per group | Yes (18 vs 18) | No (3 vs 3) | Yes if confirmed |
| Raw counts | NCBI-generated raw count matrix expected for human RNA-seq series [TO VERIFY on the Download page] | Submitter files and NCBI counts expected [TO VERIFY] | [TO VERIFY] |
| Metadata quality | Tissue type per sample, patient pairing [TO VERIFY how it is encoded] | Cell line, MOI and series are mixed in one GEO series: must select Series 1 by hand | [TO VERIFY] |
| Main confounder | All 18 patients have stage IV disease with liver metastasis; metastasis samples must be excluded | Several cell lines and infection loads in one series; small n | Subtypes, stage, possible batch [TO VERIFY] |
| Original publication | Kim SK et al. 2014, nineteen-gene risk score for colorectal cancer prognosis [checked] | Blanco-Melo D et al. Cell 2020;181(5):1036-1045 [checked] | [TO VERIFY] |
| Manageability for 1 month | High | Medium (small n, subset selection) | Medium (more metadata to understand) |

## Datasets examined and rejected

- **GSE52778** (airway smooth muscle, dexamethasone vs untreated): only 4 donors per group, paired design, and the supplementary files on the series page are FPKM and a precomputed table rather than raw counts [checked]. It is also the textbook example dataset, so it adds little to a portfolio.

## Recommendation

**A. GSE50760**, restricted to normal colon vs primary colorectal cancer. It has the largest number of replicates, a clearly defined contrast, and a known biology that makes the result easy to interpret cautiously. Its main cost is that the paired structure needs `~ patient + group` rather than a one-variable design: this is a point to confirm with the supervisor.
