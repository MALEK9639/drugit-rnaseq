**To:** Koussai Salem
**Subject:** DrugIT internship: dataset proposal for approval (GSE50760)

Dear Mr Salem,

I have finished week 1 and I am requesting approval of the dataset before starting the main analysis.

**Recommendation:** GSE50760, normal colon vs primary colorectal cancer, 18 vs 18 samples from 18 patients (liver-metastasis samples excluded). The two other candidates I compared were GSE147507 (NHBE subset, only 3 vs 3) and GSE183947 (breast tissue, more metadata to understand). GSE52778 was rejected: 4 vs 4 and no raw counts on the series page.

**Attached / linked:** repository link, RNA-seq basics note, three-dataset comparison table, completed Dataset Selection Form.

**Decision requested:** (1) Do you approve GSE50760? (2) Because the samples are paired by patient, I propose the design `~ patient + group` instead of `~ group`. Is this acceptable?

**Main risk:** all patients had stage IV disease with liver metastasis.

Best regards,
Malek Aloui
