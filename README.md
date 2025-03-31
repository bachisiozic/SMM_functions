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
# 0. Initializing the functions
source(".../SMM_functions/TumorClassification.R")
source(".../SMM_functions/SMM_score.R")
source(".../SMM_functions/crossval_cindex.R")
```

```r
# 1. Tumor Classification
class <- SMM_TnT_Cla(
  GC_input_matrix = "path/to/GenomicClassification_input_matrix_SMM_final.txt",
  cyto = "path/to/new_cyto_SMM_final.txt",
  assembly = "hg38"
)
```
```r
head(class$class)
#    sample class.new
#1 Sample_1     tumor
#2 Sample_2 non_tumor
#3 Sample_3     tumor
#4 Sample_4     tumor
#5 Sample_5     tumor
#6 Sample_6     tumor
```

```r
head(class$matrix)
#    sample CNV_chr3.gain CNV_chr5.gain CNV_chr7.gain CNV_chr9.gain CNV_chr11.gain CNV_chr15.gain CNV_chr19.gain CNV_chr21.gain CNV_chr18.gain
#1 Sample_1             1             1             0             1              1              1              1              0              0
#2 Sample_2             0             0             0             1              0              0              0              0              0
#3 Sample_3             0             0             0             0              0              0              0              0              0
#4 Sample_4             1             1             1             1              1              1              1              1              0
#5 Sample_5             1             1             0             1              0              1              1              0              0
```

```r
# 2. Genomic Score Calculation
clinic <- read.delim("path/to/clinic_SMM_bac.txt", stringsAsFactors = FALSE)
CLINIC <- clinic[complete.cases(clinic$pfs_code) & clinic$treatment == "no", ]
head(CLINIC)
#    sample pfs_code pfs_time disease_stage imwg radiology treatment cohort seq
#1 Sample_1        0       14           SMM    2         1       yes     C1 wgs
#2 Sample_2        0     1818          MGUS MGUS         1        no     C2 wes
#3 Sample_3        0     1178           SMM    2         1        no     C2 wes
#4 Sample_4        0     1316           SMM    1         1        no     C2 wes
#5 Sample_5        0      353           SMM    0         1        no     C2 wes
#6 Sample_6        0     3106           SMM    1         1        no     C2 wes

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
```

```r
head(score$score_df)
#    sample genomic_score disease_stage cohort imwg
#1 Sample_1             2           SMM     C1    2
#2 Sample_2             0           SMM     C2    1
#3 Sample_3             0           SMM     C2    0
#4 Sample_4             0           SMM     C2    1
#5 Sample_5             2          MGUS     C2 MGUS
#6 Sample_6             0           SMM     C2    1
```

```r
# 3. Cross-validation
crossval_cindex(
  data = score$analysis_data,
  output_dir = "SMM_GenomicScore/",
  k_folds = 5,
  fdr_threshold = 0.1
)
# $cindex_values
# [1] 0.6258993 0.6760204 0.7547771 0.7770701 0.6490683

# $selected_features
# sig_feat_list
#          CNV.Sig              MYC         SNV_NRAS CNV_Del_10q24.32      CNV.SNV_NF1     CNV.SNV_TET2              seq        SNV_FGFR3           APOBEC 
#                5                5                5                3                3                3                3                2                1 
#   CNV_chr15.gain   CNV_Del_2q31.1  CNV_Del_8q24.21   CNV.SNV_CREBBP   CNV.SNV_DNMT3A    CNV.SNV_NCOR1     CNV.SNV_POT1   CNV.SNV_TENT5C 
#                1                1                1                1                1                1                1                1 
```

## 📂 Output Directory

Expected files in SMM_GenomicScore/:

- `genomic_feature_stats.txt`: Genomic feature matrix (CNVs, SNVs, etc.)
- `/plots/HR_barplot.pdf, HR_with_CI_plot.pdf, crossval_cindex_boxplot.pdf`: Cytogenetic feature matrix
- `/km_plots/`: KM plots per feature
- `/cox_analysis/`: Cox model summary and forest plot

---

## 📋 Requirements

This package depends on:

- `survival`
- `survminer`
- `forestmodel`
- `caret`

These will be installed automatically if missing.

---

## Citation
Please cite our work if you use this package in your research. [Add citation here]

---

## 📬 Contact

For questions or contributions, please open an issue or contact [your email / GitHub username].






