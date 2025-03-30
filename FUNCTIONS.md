# 📘 Function Documentation

This document describes the three core functions in the **SMM Genomic Classifier & Scoring** pipeline.

---

## Function: `SMM_TnT_Cla()`

**Description**  
Classifies samples as `tumor` or `non_tumor` based on genomic and cytogenetic features, using rules informed by multiple myeloma progression biology.

```r
result <- SMM_TnT_Cla(GC_input_matrix, cyto, assembly = "hg38")`
```

**Arguments**

- `GC_input_matrix`: `data.frame` or `character`  
  A genomic feature matrix or path to a tab-delimited file containing CNVs, SNVs, etc. Must contain exactly 133 data columns.

- `cyto`: `data.frame` or `character`  
  A matrix or file containing cytogenetic information. Required columns: `sample`, `CCND3`, `CCND2`, `MAFA`, `MYC`, `hy`.

- `assembly`: `character`  
  Genome build version. Accepts `"hg38"` or `"hg19"`. Determines the set of canonical features.

**Returns**

- A named list with:
  - `class`: A data frame with tumor classification (`sample`, `class.new`).
  - `matrix`: A merged genomic matrix with input and derived features.

---

## Function: `genomic_score()`

**Description**  
Calculates a genomic score per sample using survival analysis (Cox regression and Kaplan-Meier) to identify progression-associated genomic features.

```r
score <- genomic_score(MATRIX, CLINIC, CLASS, output_dir = "results/")`
```

**Arguments**

- `MATRIX`: `data.frame` or `character`  
  Genomic matrix output from `SMM_TnT_Cla()` or a file path.

- `CLINIC`: `data.frame` or `character`  
  Clinical data frame or file path containing columns: `sample`, `pfs_time`, `pfs_code`, `disease_stage`, `imwg`.

- `CLASS`: `data.frame` or `character`  
  Tumor classification (from `SMM_TnT_Cla()`), optional but recommended.

- `output_dir`: `character`  
  Output directory where all plots and tables will be saved.

- `filter_value`: `numeric`  
  Threshold for selecting significant features (default: `0.1` for FDR).

- `selection_by`: `character`  
  Filtering method: `"all"`, `"stage"`, or `"cohort"`.

- `selection_value`: `character` or `vector`  
  Value(s) corresponding to the `selection_by` method.

- `filter_by`: `character`  
  Filtering metric: `"FDR"` or `"p_value"`.

- `min_event_count`: `numeric`  
  Minimum number of samples in which a feature must be present (default: `2`).

- `use_known_features`: `logical`  
  Use pre-defined features instead of testing all (default: `FALSE`).

**Returns**

- A named list with:
  - `score_df`: A data frame with per-sample genomic scores and metadata.
  - `analysis_data`: Merged dataset used for statistical analysis.

**Output Files (in `output_dir`)**

- `genomic_feature_stats.txt`: Stats for all tested features (HR, p-value, FDR, counts).
- `genomic_score_per_sample.txt`: Final scores by sample.
- `/plots/`: Barplots and forest plots of hazard ratios.
- `/km_plots/`: Kaplan-Meier survival curves for significant features.
- `/cox_analysis/`: Multivariate Cox model results and forest plot.

---

## Function: `crossval_cindex()`

**Description**  
Performs `k`-fold cross-validation using feature selection and Cox models. Reports predictive performance via the concordance index (`c-index`).

 
```r
crossval_cindex(data = score$analysis_data, output_dir = "results/", k_folds = 5, fdr_threshold = 0.1)`
```

**Arguments**

- `data`: `data.frame`  
  Data object from `genomic_score()$analysis_data`.

- `output_dir`: `character`  
  Output directory where results (plots) will be saved.

- `k_folds`: `numeric`  
  Number of cross-validation folds (default: `5`).

- `fdr_threshold`: `numeric`  
  FDR threshold used to select features in each training set (default: `0.1`).

**Returns**

- A named list with:
  - `cindex_values`: Numeric vector of concordance index values for each fold.
  - `selected_features`: Ranked table of selected features across all folds.

**Output Files**

- `/plots/crossval_cindex_boxplot.pdf`: Boxplot of c-index values across folds.

---

## 📦 Dependencies

The following R packages are required and will be loaded automatically by the functions:

- `survival`
- `survminer`
- `forestmodel`
- `caret`

These will be installed automatically if not already present in your environment.

---

## 📬 Contact

For questions or contributions, please contact the package maintainer or open an issue in the GitHub repository.
