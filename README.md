# 🔬 SMM Genomic Classifier

This R-based toolset includes two key functions to assist defining malignant transformation begins with the analysis of a matrix comprising 136 genomic features (plus MYC and RAS pathway mutations), including aneuploidies, mutations, and translocations, extracted from patient samples.

---

## 📦 Included Functions

1. **`genomic_score()`** – Compute the genomic drivers associated with progression.
2. **`SMM_TnT_Cla()`** – Classify tumor vs. non-tumor samples from genomic features.

Help in `FUNCTIONS.md`

---

## 🧪 Input Files

- `GC_input_matrix`: Genomic feature matrix (CNVs, SNVs, etc.) [https://github.com/UM-Myeloma-Genomics/GCP_MM]
- `cyto`: Cytogenetic feature matrix
- `CLINIC`: Clinical data with `pfs_time`, `pfs_code`, `stage`, etc.
- `MATRIX`: Merged genomic matrix output from classification

---

## 🔁 Pipeline Steps

```r
# 0. Initializing the functions
source(".../SMM_functions/1.TumorClassification.R")
source(".../SMM_functions/2.SMM_score.R")
```

```r
# 1. Prognostic drivers
training.mtx=read.delim("~/Desktop/BACHI/Project1/MAY2025/0.Final/final_matrix_training.txt",
                        stringsAsFactors = F)
training.clin=read.delim("~/Desktop/BACHI/Project1/MAY2025/0.Final/final_clinic_training.txt",
                         stringsAsFactors = F)
#    sample pfs_code pfs_time disease_stage imwg radiology treatment cohort seq
#1 Sample_1        0       14           SMM    2         1       yes     C1 wgs
#2 Sample_2        0     1818          MGUS MGUS         1        no     C2 wes
#3 Sample_3        0     1178           SMM    2         1        no     C2 wes
#4 Sample_4        0     1316           SMM    1         1        no     C2 wes
#5 Sample_5        0      353           SMM    0         1        no     C2 wes
#6 Sample_6        0     3106           SMM    1         1        no     C2 wes

training.mtx$RAS=rowSums(training.mtx[,c(68,71,76,78,80,52)]) # BRAF, KRAS, NRAS, FGFR3, NF1, PTPN11
training.mtx$RAS=ifelse(training.mtx$RAS >1,1,training.mtx$RAS)

score=genomic_score(training.mtx, 
                    training.clin, 
                    selection_by = "cohort", 
                    selection_value = c("MSKCC","bolli","Sanger","Moffitt"), 
                    filter_by = "p.value",
                    filter_value = 0.05,
                    min_event_count=3,
                    output_dir="~/Desktop/BACHI/Project1/MAY2025/SMM_GenomicScore/")

score$sig_features
#[1] "CNV_chr3.gain"    "CNV_chr5.gain"    "CNV_chr9.gain"    "CNV_chr11.gain"   "CNV_chr15.gain"   "CNV_chr19.gain"   "CNV_chr21.gain"   "CNV_Del_10q24.32"
#[9] "CNV_Del_2q31.1"   "CNV_Del_6q26"     "CNV_Del_8q24.21"  "CNV.SNV_ARID2"    "CNV.SNV_CREBBP"   "CNV.SNV_CYLD"     "CNV.SNV_DNMT3A"   "CNV.SNV_TENT5C"  
#[17] "CNV.SNV_NCOR1"    "CNV.SNV_NF1"      "CNV.SNV_POT1"     "CNV.SNV_PRDM1"    "CNV.SNV_TET2"     "SNV_FGFR3"        "SNV_NRAS"         "CNV.Sig"         
#[25] "APOBEC"           "t_MMSET"          "MYC"              "hy"               "RAS"   
```

```r
# 2. Tumor Classification
class=SMM_TnT_Cla(GC_input_matrix = "~/Desktop/BACHI/Project1/MAY2025/0.Final/final_genomic_matrix_374pts.txt",
                  cyto = "~/Desktop/BACHI/Project1/MAY2025/0.Final/final_cytogenetic_374pts.txt",
                  sig.features = c(score$sig_features))
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

For questions or contributions, please open an issue or contact Bachisio Ziccheddu (bxz262@miami.edu) & Francesco Maura (mauraf@mskcc.org).






