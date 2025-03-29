crossval_cindex <- function(data, output_dir = "genomic_score_output",k_folds=5, fdr_threshold = 0.1) {
  required_packages <- c("caret", "survival", "survminer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  set.seed(123)
  
  folds <- createFolds(data$pfs_code, k = k_folds, list = TRUE, returnTrain = FALSE)
  cindex_values <- c()
  sig_feat_list <- c()
  
  for (index in folds) {
    train_data <- data[-index, ]
    test_data <- data[index, ]
    
    tumor_fold_train <- list()
    
    for (i in 1:ncol(train_data)) {
      feature_name <- colnames(train_data)[i]
      if (feature_name %in% c("sample", "pfs_time", "pfs_code", "disease_stage", "cohort", "imwg", "class")) next
      
      train_data$test <- train_data[[feature_name]]
      if (length(unique(train_data$test)) < 2 || sum(train_data$test >= 1, na.rm=TRUE) < 3) next
      
      if (length(unique(train_data$test)) == 3 && all(unique(train_data$test) %in% c(0, 1, 2))) {
        train_data$test[train_data$test > 1] <- 1
      }
      
      km_result <- survminer::surv_fit(survival::Surv(pfs_time, pfs_code) ~ test, data = train_data)
      km_p_value <- as.numeric(gsub("p = ", "", survminer::surv_pvalue(km_result)[, 2]))
      
      fit <- survival::coxph(Surv(pfs_time, pfs_code) ~ test, data = train_data)
      hr <- exp(coef(fit))
      ci <- exp(confint(fit))
      
      tumor_fold_train[[feature_name]] <- c(feature_name, km_p_value, hr, ci[1], ci[2], sum(train_data$test == 1))
    }
    
    tumor_fold_train_df <- do.call(rbind.data.frame, tumor_fold_train)
    if (nrow(tumor_fold_train_df) == 0) next
    
    colnames(tumor_fold_train_df) <- c("driver", "p_value", "HR", "HR_low", "HR_high", "present")
    tumor_fold_train_df$p_value <- as.numeric(tumor_fold_train_df$p_value)
    tumor_fold_train_df$FDR <- p.adjust(tumor_fold_train_df$p_value, method = "fdr")
    
    selected_features <- tumor_fold_train_df$driver[tumor_fold_train_df$FDR < fdr_threshold]
    sig_feat_list <- c(sig_feat_list, selected_features)
    
    if (length(selected_features) > 0) {
      cox_model <- coxph(Surv(pfs_time, pfs_code) ~ ., data = train_data[, c("pfs_time", "pfs_code", "imwg", selected_features), drop=FALSE])
      test_data$predicted_risk <- predict(cox_model, newdata = test_data)
      validation_model <- coxph(Surv(pfs_time, pfs_code) ~ predicted_risk, data = test_data)
      c_index <- concordance(validation_model)$concordance
      cindex_values <- c(cindex_values, c_index)
    }
  }
  
  # C-index boxplot
  pdf(file.path(output_dir, "plots", "crossval_cindex_boxplot.pdf"), width = 4, height = 6)
  boxplot(x = na.omit(as.numeric(cindex_values)), 
          las = 2, ylab = "c-index", xlab = "", xaxt = "n", bty = "n", pch = 16, xlim = c(0.8, 1.2))
  par(new = TRUE)
  plot(y = na.omit(as.numeric(cindex_values)), x = rep(1, length(cindex_values)), las = 2, ylab = "c-index",
       xlab = "", xaxt = "n", bty = "n", pch = 21, xlim = c(0.8, 1.2), bg = "coral", cex = 2)
  segments(y0 = min(na.omit(as.numeric(cindex_values))),
           y1 = max(na.omit(as.numeric(cindex_values))), 
           x0 = rep(1, length(cindex_values)),
           x1 = rep(1, length(cindex_values)))
  segments(y0 = median(na.omit(as.numeric(cindex_values))),
           y1 = median(na.omit(as.numeric(cindex_values))), 
           x0 = 0.9,
           x1 = 1.1)
  dev.off()
  
  return(list(cindex_values = cindex_values, selected_features = sort(table(sig_feat_list), decreasing = TRUE)))
}
