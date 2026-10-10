# RNA-seq basics (one-page note)

*Draft: rewrite it in your own words before sending. The assignment asks you to show that you understand the concepts.*

**Gene and transcript.** A gene is a DNA region that can be transcribed. One gene can produce several transcripts (isoforms) through alternative splicing. Gene-level analysis sums the signal of all isoforms of a gene.

**Read.** A short DNA/cDNA sequence (typically 50-150 bases) produced by the sequencer. Reads are aligned to a reference genome (or pseudo-aligned to transcripts) to find the gene they come from.

**Count matrix.** A table with genes as rows and samples as columns; each cell is the number of reads assigned to that gene in that sample. Raw counts are integers.

**Sample metadata.** A table with one row per sample: sample ID, group (e.g. normal / tumor), and any other variable (patient, batch, sex). The sample IDs must match the count-matrix column names exactly.

**Biological replicate.** An independent biological unit in the same condition (different patient, different culture). Replicates let us estimate biological variability, which is what makes a statistical test possible. Sequencing the same library twice is a technical replicate and does not replace biological replicates.

**Raw counts, normalized counts, TPM, FPKM.**
- *Raw counts*: reads per gene, depend on sequencing depth. Required by DESeq2, which models the counting noise itself.
- *Normalized counts*: raw counts divided by sample-specific size factors (DESeq2 median-of-ratios) so samples are comparable.
- *TPM / FPKM*: expression scaled by gene length and library size, useful to compare genes within a sample. They are not valid input for DESeq2.

**Differential expression.** A statistical test, gene by gene, of whether the mean expression differs between two groups more than expected from biological variability. The output is a log2 fold change (effect size), a p-value, and an adjusted p-value that corrects for testing thousands of genes at once (Benjamini-Hochberg, false discovery rate).

**Why a paired design matters here.** When the same patient provides a normal and a tumour sample, the patient is a known source of variation; including it in the model (`~ patient + group`) compares each tumour with its own normal tissue.

**Wording reminder.** Say "differentially expressed in this dataset", not "causes"; "candidate expression signal requiring validation", not "biomarker discovered".
