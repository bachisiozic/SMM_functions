# 🔬 SMM Genomic Classifier & Scoring Pipeline

This R-based toolset includes three key functions to assist in classifying tumor samples, scoring genomic risk in smoldering multiple myeloma (SMM), and validating prognostic performance via cross-validation.

---

## 📦 Included Functions

1. **`SMM_TnT_Cla()`** – Classify tumor vs. non-tumor samples from genomic features.
2. **`genomic_score()`** – Calculate genomic scores predictive of progression.
3. **`crossval_cindex()`** – Perform cross-validation to assess prognostic accuracy.

---

## 🧪 Input Files

- `GC_input_matrix`: Genomic feature matrix (CNVs, SNVs, etc.)
- `cyto`: Cytogenetic feature matrix
- `CLINIC`: Clinical data with `pfs_time`, `pfs_code`, `stage`, etc.
- `CLASS`: Output from classification script
- `MATRIX`: Merged genomic matrix output from classification

---

## 🔁 Pipeline Steps

```r
# 1. Tumor Classification
class <- SMM_TnT_Cla(
  GC_input_matrix = "path/to/GenomicClassification_input_matrix_SMM_final.txt",
  cyto = "path/to/new_cyto_SMM_final.txt",
  assembly = "hg38"
)

# 2. Genomic Score Calculation
clinic <- read.delim("path/to/clinic_SMM_bac.txt", stringsAsFactors = FALSE)
CLINIC <- clinic[complete.cases(clinic$pfs_code) & clinic$treatment == "no", ]
MATRIX <- class$matrix[class$matrix$sample %in% CLINIC$sample, ]
CLASS <- class$class[class$class$sample %in% CLINIC$sample, ]

score <- genomic_score(
  MATRIX, CLINIC, CLASS,
  selection_by = "cohort",
  selection_value = c("MSKCC", "bolli", "Sanger", "Moffitt"),
  filter_by = "FDR",
  filter_value = 0.1,
  min_event_count = 3,
  output_dir = "SMM_GenomicScore/"
)

# 3. Cross-validation
crossval_cindex(
  data = score$analysis_data,
  output_dir = "SMM_GenomicScore/",
  k_folds = 5,
  fdr_threshold = 0.1
)
```
