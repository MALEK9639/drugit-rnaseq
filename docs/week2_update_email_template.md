**Subject:** DrugIT internship: week 2 update (data preparation and QC)

**Completed this week**
- Downloaded raw counts and metadata for GSE50760 (log: `metadata/retrieval_log.csv`)
- Cleaned metadata and counts (log: `metadata/cleaning_log.txt`)
- Quality control: library size, filtering, sample distances, PCA

**Repository:** <link>

**Main figures, one interpretation each** (fill in from your real results)
- `figures/library_size.png`: <one sentence, e.g. library sizes range from X M to Y M, no sample below half of the median>
- `figures/pca_plot.png`: <PC1 explains X%; do the two groups separate? Do the same patients cluster?>
- `figures/sample_distance_heatmap.png`: <do samples cluster by group?>

**Filtering rule:** count >= 10 in at least <n> samples; <before> genes before, <after> after.
**Missing values:** <count and action>.
**Outliers:** <flagged samples from `results/qc_outlier_flags.csv`, written decision, no removal without your agreement>.

**Main difficulty:** <one item>
**Question requiring your input:** <only if needed, using the technical-question format>
**Plan for next week:** DESeq2 with the approved design, MA and volcano plots, top-gene heatmap, pathway enrichment.
