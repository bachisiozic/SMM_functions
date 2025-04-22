genomic_score=function(MATRIX, CLINIC, CLASS = NULL, output_dir = "genomic_score_output", filter_value = 0.1,
                          selection_by = "all", selection_value = NULL,
                          filter_by = "FDR", min_event_count = 2, use_known_features = FALSE) {
  # Required packages
  required_packages=c("survival", "forestmodel", "survminer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Create output subdirectories if they don't exist
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "km_plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "cox_analysis"), showWarnings = FALSE, recursive = TRUE)
  
  if (class(MATRIX)=="character") {
    matrix=read.delim(MATRIX, stringsAsFactors=FALSE)
  } else {
    matrix=MATRIX
  }
  
  if (class(CLINIC)=="character") {
    clinic=read.delim(CLINIC, stringsAsFactors=FALSE)
  } else {
    clinic=CLINIC
  }
  
  required_cols=c("sample", "pfs_code", "pfs_time", "disease_stage", "imwg")
  if (!all(required_cols %in% colnames(clinic))) {
    stop("Missing required clinical columns: sample, pfs_code, pfs_time, disease_stage, imwg")
  }
  
  ## Arranging MATRIX
  data=MATRIX
  copynumbers=colnames(data)[grep("CNV_",colnames(data))]
  tumorsuppressorgenes=colnames(data)[grep("CNV.SNV_",colnames(data))]
  oncogenesothers=c(colnames(data)[grep("^SNV_|^SNV.",colnames(data))])
  genomics=data[,c(colnames(data)[grep("sample|sampleID|Sample|SampleID|SAMPLE|SAMPLE.ID",colnames(data))],
                   copynumbers,tumorsuppressorgenes,oncogenesothers,
                   "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy")]
  
  ## Merge
  all_data=merge(clinic, matrix, by="sample")
  
  # Merge CLASS if provided
  if (!is.null(CLASS)) {
    if (class(CLASS)=="character") {
      class_df=read.delim(CLASS, stringsAsFactors=FALSE)
    } else {
      class_df=CLASS
    }
  }
  
  if (!is.null(CLASS)) {
    if (!"sample" %in% colnames(class_df)) stop("CLASS must contain a 'sample' column")
    all_data=merge(all_data, class_df, by = "sample", all.x = TRUE)
    colnames(all_data)[ncol(all_data)]="class"
  }
  
  
  # Always apply tumor filtering first if CLASS is present
  if (!is.null(CLASS) && "class" %in% colnames(all_data)) {
    removed_classes=table(all_data$class)
    all_data=all_data[all_data$class == "tumor", ]
    cat("Tumor filtering applied.")
    cat("Samples removed by class filter:")
    print(removed_classes[names(removed_classes) != "tumor"])
    cat("Remaining tumor samples:", nrow(all_data), "")
  }
  #}  # Close tumor filtering and then CLASS
  #}
  #  ]
  #}  # End of CLASS merge and tumor filtering
  
  # Determine selection set
  cat("=== SELECTION FILTERING REPORT ===\n")
  initial_n=nrow(all_data)
  if (selection_by == "stage") {
    data=all_data[all_data$disease_stage %in% selection_value & !is.na(all_data$pfs_code), ]
    cat("Selection by disease stage:", selection_value,"\n",sep="")
    cat("Samples removed by stage filter:", initial_n - nrow(data),"\n" ,sep="")
  } else if (selection_by == "cohort") {
    data=all_data[all_data$cohort %in% selection_value & !is.na(all_data$pfs_code), ]
    cat("Selection by cohort:", selection_value,"\n",sep=" ")
    cat("Samples removed by cohort filter:", initial_n - nrow(data),"\n",sep="")
  } else {
    data=all_data[all_data$disease_stage %in% c("SMM", "MGUS") & !is.na(all_data$pfs_code), ]
    cat("Default selection (SMM and MGUS).")
    cat("Samples removed by default stage filter:", initial_n - nrow(data), "")
  }
  
  
  if (use_known_features) {
    sig_features=c("CNV.SNV_TET2", "CNV.SNV_NF1", "CNV_Del_10q24.32", "MYC", "CNV.Sig", "SNV_NRAS")
    sig_features=sig_features[sig_features %in% colnames(all_data)]
    stat_df=data.frame()
    for (feature in sig_features) {
      x=data
      x$feature=x[[feature]]
      if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
        x$feature[x$feature > 1]=1
      }
      if (length(unique(x$feature)) < 2) next
      fit=coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
      hr=exp(coef(fit))
      ci=exp(confint(fit))
      result_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
      km_p_value=as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
      pval=km_p_value
      present=sum(x$feature %in% c(1, 2), na.rm=TRUE)
      absent=sum(x$feature == 0, na.rm=TRUE)
      present_MGUS=sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
      present_SMM=sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
      stat_df=rbind(stat_df, data.frame(
        feature = feature,
        p_value = km_p_value,
        HR = hr,
        HR_low = ci[1],
        HR_high = ci[2],
        present = present,
        absent = absent,
        present_MGUS = present_MGUS,
        present_SMM = present_SMM,
        FDR = NA
      ))
    }
    stat_df$FDR=p.adjust(stat_df$p_value, method="fdr")
    write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
  } else {
    stat_df=data.frame()
    #print(colnames(data)[!colnames(data) %in% c("sample", "pfs_time", "pfs_code", "treatment","disease_stage", "cohort", "imwg", "class", "seq", "radiology")])
    for (feature in colnames(data)[!colnames(data) %in% c("sample", "pfs_time", "pfs_code", "treatment","disease_stage", "cohort", "imwg", "class", "seq", "radiology")]) {
      x=data
      x$feature=x[[feature]]
      if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
        x$feature[x$feature > 1]=1
      }
      if (length(unique(x$feature)) < 2) next
      fit=coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
      hr=exp(coef(fit))
      ci=exp(confint(fit))
      result_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
      km_p_value=as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
      pval=km_p_value
      present=sum(x$feature %in% c(1, 2), na.rm=TRUE)
      absent=sum(x$feature == 0, na.rm=TRUE)
      present_MGUS=sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
      present_SMM=sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
      stat_df=rbind(stat_df, data.frame(
        feature = feature,
        p_value = pval,
        HR = hr,
        HR_low = ci[1],
        HR_high = ci[2],
        present = present,
        absent = absent,
        present_MGUS = present_MGUS,
        present_SMM = present_SMM,
        FDR = NA
      ))
    }
    stat_df$FDR=as.numeric(as.character(stat_df$FDR))
    stat_df$p_value=as.numeric(as.character(stat_df$p_value))
    
    stat_df$NOTE=ifelse(stat_df$present < min_event_count,"Not enough events",".")
    stat_df1=stat_df[stat_df$NOTE != "Not enough events",]
    stat_df2=stat_df[stat_df$NOTE == "Not enough events",]
    
    stat_df1$FDR=p.adjust(as.numeric(as.character(stat_df1$p_value)), method = "fdr")
    stat_df=rbind(stat_df1,stat_df2)
    
    stat_df=stat_df[match(c(copynumbers,tumorsuppressorgenes,oncogenesothers,
                            "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy"),
                          stat_df$feature),]
    stat_df=stat_df[complete.cases(stat_df$feature),]
    write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
    
    if (filter_by == "FDR") {
      sig_features1=na.omit(stat_df$feature[stat_df$FDR < filter_value & stat_df$present >= min_event_count])
    } else {
      sig_features1=na.omit(stat_df$feature[stat_df$p_value < filter_value & stat_df$present >= min_event_count])
    }
    print(na.omit(sig_features1))
  }
  
  if (length(sig_features1) > 0) {
    all_data$genomic_score=rowSums(all_data[, sig_features1, drop=FALSE], na.rm=TRUE)
  } else {
    warning("None of the selected features are present in the input data. Setting all genomic scores to 0.")
    all_data$genomic_score=0
  }
  all_data$genomic_score[all_data$genomic_score > 0]=2
  
  
  score_df=data.frame(sample = all_data$sample,
                         genomic_score = all_data$genomic_score,
                         disease_stage = all_data$disease_stage,
                         cohort = all_data$cohort,
                         imwg = all_data$imwg)
  write.table(score_df, file.path(output_dir, "genomic_score_per_sample.txt"), sep="	", row.names=FALSE, quote=FALSE)
  
  # Barplot of HR
  # pdf(file.path(output_dir, "plots", "HR_barplot.pdf"), onefile = TRUE)
  # hr_data=stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
  # if (nrow(hr_data) > 0) {
  #   barplot(hr_data$HR,
  #           names.arg = hr_data$feature,
  #           las = 2, col = "coral", ylab = "Hazard Ratio")
  # } else {
  #   plot.new()
  #   text(0.5, 0.5, "No finite HR values with >=3 events to plot")
  # }
  # dev.off()
  
  # # Plot HR with CI error bars
  # pdf(file.path(output_dir, "plots", "HR_with_CI_plot.pdf"), onefile = TRUE)
  # stat_df_sig=stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
  # if (nrow(stat_df_sig) > 0 && any(is.finite(stat_df_sig$HR_high))) {
  #   par(mar = c(10, 5, 4, 2))
  #   ylim_max=max(stat_df_sig$HR_high, na.rm = TRUE)
  #   if (is.finite(ylim_max)) {
  #     x=barplot(stat_df_sig$HR, names.arg = stat_df_sig$feature, las = 2,
  #                  col = "skyblue", ylab = "Hazard Ratio", ylim = c(0, ylim_max * 1.1))
  #     segments(x0 = x, x1 = x, y0 = stat_df_sig$HR_low, y1 = stat_df_sig$HR_high)
  #   } else {
  #     plot.new()
  #     text(0.5, 0.5, "No finite HR CI to plot")
  #   }
  # } else {
  #   plot.new()
  #   text(0.5, 0.5, "No significant features with valid HR and CI")
  # }
  # dev.off()
  
  # KM plots for significant features
  for (feature in sig_features1) {
    x=data
    x$feature=x[[feature]]
    if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
      x$feature[x$feature > 1]=1
    }
    if (length(unique(x$feature)) < 2) next
    fit_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
    pdf(file.path(output_dir, "km_plots", paste0("KM_", gsub("/", "_", feature), ".pdf")), onefile = FALSE)
    print(ggsurvplot(fit_km, data = x, pval = TRUE, risk.table = TRUE, xlab = "Time", ylab = feature))
    dev.off()
    
    
  }
  
  # Multivariate Cox model using all samples
  if (length(sig_features1) >= 1 && any(!is.na(all_data$pfs_code))) {
    full_data=all_data[!is.na(all_data$pfs_code), ]
    if (nrow(full_data) > 0) {
      cox_all=coxph(Surv(pfs_time, pfs_code) ~ genomic_score + imwg, data = full_data, robust = TRUE)
      summary_all=as.data.frame(summary(cox_all)$coefficients)
      summary_all$variable=rownames(summary_all)
      write.table(summary_all, file.path(output_dir, "cox_analysis", "cox_summary_all_samples.txt"), sep="	", row.names=FALSE, quote=FALSE)
      
      pdf(file.path(output_dir, "cox_analysis", "forest_plot_all_samples.pdf"))
      print(forest_model(cox_all))
      dev.off()
    }
  }
  
  result=list(
    score_df = score_df,
    analysis_data = data,
    sig_features = sig_features1
  )
  return(result)
}



# genomic_score=function(MATRIX, CLINIC, CLASS = NULL, output_dir = "genomic_score_output", filter_value = 0.1,
#                           selection_by = "all", selection_value = NULL,
#                           filter_by = "FDR", min_event_count = 2, use_known_features = FALSE) {
#   # Required packages
#   required_packages=c("survival", "forestmodel", "survminer")
#   for (pkg in required_packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       install.packages(pkg)
#     }
#     library(pkg, character.only = TRUE)
#   }
#   
#   # Create output subdirectories if they don't exist
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
#   dir.create(file.path(output_dir, "km_plots"), showWarnings = FALSE, recursive = TRUE)
#   dir.create(file.path(output_dir, "cox_analysis"), showWarnings = FALSE, recursive = TRUE)
#   
#   if (class(MATRIX)=="character") {
#     matrix=read.delim(MATRIX, stringsAsFactors=FALSE)
#   } else {
#     matrix=MATRIX
#   }
#   
#   if (class(CLINIC)=="character") {
#     clinic=read.delim(CLINIC, stringsAsFactors=FALSE)
#   } else {
#     clinic=CLINIC
#   }
#   
#   required_cols=c("sample", "pfs_code", "pfs_time", "disease_stage", "imwg")
#   if (!all(required_cols %in% colnames(clinic))) {
#     stop("Missing required clinical columns: sample, pfs_code, pfs_time, disease_stage, imwg")
#   }
#   
#   ## Merge
#   all_data=merge(clinic, matrix, by="sample")
#   
#   # Merge CLASS if provided
#   if (!is.null(CLASS)) {
#     if (class(CLASS)=="character") {
#       class_df=read.delim(CLASS, stringsAsFactors=FALSE)
#     } else {
#       class_df=CLASS
#     }
#   }
#   
#   if (!is.null(CLASS)) {
#     if (!"sample" %in% colnames(class_df)) stop("CLASS must contain a 'sample' column")
#     all_data=merge(all_data, class_df, by = "sample", all.x = TRUE)
#   }
#   
#   colnames(all_data)[ncol(all_data)]="class"
#   
#   # Always apply tumor filtering first if CLASS is present
#   if (!is.null(CLASS) && "class" %in% colnames(all_data)) {
#     removed_classes=table(all_data$class)
#     all_data=all_data[all_data$class == "tumor", ]
#     cat("Tumor filtering applied.
# ")
#     cat("Samples removed by class filter:
# ")
#     print(removed_classes[names(removed_classes) != "tumor"])
#     cat("Remaining tumor samples:", nrow(all_data), "
# ")
#   }
#   #}  # Close tumor filtering and then CLASS
#   #}
#   #  ]
#   #}  # End of CLASS merge and tumor filtering
#   
#   # Determine selection set
#   cat("
# === SELECTION FILTERING REPORT ===
# ")
#   initial_n=nrow(all_data)
#   if (selection_by == "stage") {
#     data=all_data[all_data$disease_stage %in% selection_value & !is.na(all_data$pfs_code), ]
#     cat("Selection by disease stage:", selection_value, "
# ")
#     cat("Samples removed by stage filter:", initial_n - nrow(data), "
# ")
#   } else if (selection_by == "cohort") {
#     data=all_data[all_data$cohort %in% selection_value & !is.na(all_data$pfs_code), ]
#     cat("Selection by cohort:", selection_value, "
# ")
#     cat("Samples removed by cohort filter:", initial_n - nrow(data), "
# ")
#   } else {
#     data=all_data[all_data$disease_stage %in% c("SMM", "MGUS") & !is.na(all_data$pfs_code), ]
#     cat("Default selection (SMM and MGUS).
# ")
#     cat("Samples removed by default stage filter:", initial_n - nrow(data), "
# ")
#   }
#   
#   
#   
#   
#   if (use_known_features) {
#     sig_features=c("CNV.SNV_TET2", "CNV.SNV_NF1", "CNV_Del_10q24.32", "MYC", "CNV.Sig", "SNV_NRAS")
#     sig_features=sig_features[sig_features %in% colnames(all_data)]
#     stat_df=data.frame()
#     for (feature in sig_features) {
#       x=data
#       x$feature=x[[feature]]
#       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#         x$feature[x$feature > 1]=1
#       }
#       if (length(unique(x$feature)) < 2) next
#       fit=coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       hr=exp(coef(fit))
#       ci=exp(confint(fit))
#       result_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       km_p_value=as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
#       pval=km_p_value
#       present=sum(x$feature %in% c(1, 2), na.rm=TRUE)
#       absent=sum(x$feature == 0, na.rm=TRUE)
#       present_MGUS=sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
#       present_SMM=sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
#       stat_df=rbind(stat_df, data.frame(
#         feature = feature,
#         p_value = km_p_value,
#         HR = hr,
#         HR_low = ci[1],
#         HR_high = ci[2],
#         present = present,
#         absent = absent,
#         present_MGUS = present_MGUS,
#         present_SMM = present_SMM,
#         FDR = NA
#       ))
#     }
#     stat_df$FDR=p.adjust(stat_df$p_value, method="fdr")
#     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
#   } else {
#     stat_df=data.frame()
#     for (feature in colnames(data)[!colnames(data) %in% c("sample", "pfs_time", "pfs_code", "disease_stage", "cohort", "imwg", "class", "seq", "hy")]) {
#       x=data
#       x$feature=x[[feature]]
#       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#         x$feature[x$feature > 1]=1
#       }
#       if (length(unique(x$feature)) < 2) next
#       fit=coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       hr=exp(coef(fit))
#       ci=exp(confint(fit))
#       result_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       km_p_value=as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
#       pval=km_p_value
#       present=sum(x$feature %in% c(1, 2), na.rm=TRUE)
#       absent=sum(x$feature == 0, na.rm=TRUE)
#       present_MGUS=sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
#       present_SMM=sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
#       stat_df=rbind(stat_df, data.frame(
#         feature = feature,
#         p_value = pval,
#         HR = hr,
#         HR_low = ci[1],
#         HR_high = ci[2],
#         present = present,
#         absent = absent,
#         present_MGUS = present_MGUS,
#         present_SMM = present_SMM,
#         FDR = NA
#       ))
#     }
#     stat_df$FDR=p.adjust(stat_df$p_value, method = "fdr")
#     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
#     
#     if (filter_by == "FDR") {
#       sig_features=stat_df$feature[stat_df$FDR < filter_value & stat_df$present >= min_event_count]
#     } else {
#       sig_features=stat_df$feature[stat_df$p_value < 0.05 & stat_df$present >= min_event_count]
#     }
#     sig_features=sig_features[sig_features %in% colnames(all_data)]
#   }
#   
#   if (length(sig_features) > 0) {
#     all_data$genomic_score=rowSums(all_data[, sig_features, drop=FALSE], na.rm=TRUE)
#   } else {
#     warning("None of the selected features are present in the input data. Setting all genomic scores to 0.")
#     all_data$genomic_score=0
#   }
#   all_data$genomic_score[all_data$genomic_score > 0]=2
#   
#   
#   score_df=data.frame(sample = all_data$sample,
#                          genomic_score = all_data$genomic_score,
#                          disease_stage = all_data$disease_stage,
#                          cohort = all_data$cohort,
#                          imwg = all_data$imwg)
#   write.table(score_df, file.path(output_dir, "genomic_score_per_sample.txt"), sep="	", row.names=FALSE, quote=FALSE)
#   
#   # Barplot of HR
#   pdf(file.path(output_dir, "plots", "HR_barplot.pdf"), onefile = TRUE)
#   hr_data=stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
#   if (nrow(hr_data) > 0) {
#     barplot(hr_data$HR,
#             names.arg = hr_data$feature,
#             las = 2, col = "coral", ylab = "Hazard Ratio")
#   } else {
#     plot.new()
#     text(0.5, 0.5, "No finite HR values with >=3 events to plot")
#   }
#   dev.off()
#   
#   # Plot HR with CI error bars
#   pdf(file.path(output_dir, "plots", "HR_with_CI_plot.pdf"), onefile = TRUE)
#   stat_df_sig=stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
#   if (nrow(stat_df_sig) > 0 && any(is.finite(stat_df_sig$HR_high))) {
#     par(mar = c(10, 5, 4, 2))
#     ylim_max=max(stat_df_sig$HR_high, na.rm = TRUE)
#     if (is.finite(ylim_max)) {
#       x=barplot(stat_df_sig$HR, names.arg = stat_df_sig$feature, las = 2,
#                    col = "skyblue", ylab = "Hazard Ratio", ylim = c(0, ylim_max * 1.1))
#       segments(x0 = x, x1 = x, y0 = stat_df_sig$HR_low, y1 = stat_df_sig$HR_high)
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No finite HR CI to plot")
#     }
#   } else {
#     plot.new()
#     text(0.5, 0.5, "No significant features with valid HR and CI")
#   }
#   dev.off()
#   
#   # KM plots for significant features
#   for (feature in sig_features) {
#     x=data
#     x$feature=x[[feature]]
#     if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#       x$feature[x$feature > 1]=1
#     }
#     if (length(unique(x$feature)) < 2) next
#     fit_km=surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#     pdf(file.path(output_dir, "km_plots", paste0("KM_", gsub("/", "_", feature), ".pdf")), onefile = FALSE)
#     print(ggsurvplot(fit_km, data = x, pval = TRUE, risk.table = TRUE, xlab = "Time", ylab = feature))
#     dev.off()
#     
#     
#   }
#   
#   # Multivariate Cox model using all samples
#   if (length(sig_features) >= 1 && any(!is.na(all_data$pfs_code))) {
#     full_data=all_data[!is.na(all_data$pfs_code), ]
#     if (nrow(full_data) > 0) {
#       cox_all=coxph(Surv(pfs_time, pfs_code) ~ genomic_score + imwg, data = full_data, robust = TRUE)
#       summary_all=as.data.frame(summary(cox_all)$coefficients)
#       summary_all$variable=rownames(summary_all)
#       write.table(summary_all, file.path(output_dir, "cox_analysis", "cox_summary_all_samples.txt"), sep="	", row.names=FALSE, quote=FALSE)
#       
#       pdf(file.path(output_dir, "cox_analysis", "forest_plot_all_samples.pdf"))
#       print(forest_model(cox_all))
#       dev.off()
#     }
#   }
#   
#   return(score_df)
# }





















# #### . . . LAST


# genomic_score <- function(MATRIX, CLINIC, CLASS = NULL, output_dir = "genomic_score_output", filter_value = 0.1,
#                           selection_by = "all", selection_value = NULL,
#                           filter_by = "FDR", min_event_count = 2, use_known_features = FALSE) {
#   # Required packages
#   required_packages <- c("survival", "forestmodel", "survminer")
#   for (pkg in required_packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       install.packages(pkg)
#     }
#     library(pkg, character.only = TRUE)
#   }
  
#   # Create output subdirectories if they don't exist
#   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
#   dir.create(file.path(output_dir, "km_plots"), showWarnings = FALSE, recursive = TRUE)
#   dir.create(file.path(output_dir, "cox_analysis"), showWarnings = FALSE, recursive = TRUE)
  
#   if (class(MATRIX)=="character") {
#     matrix <- read.delim(MATRIX, stringsAsFactors=FALSE)
#   } else {
#     matrix <- MATRIX
#   }
  
#   if (class(CLINIC)=="character") {
#     clinic <- read.delim(CLINIC, stringsAsFactors=FALSE)
#   } else {
#     clinic <- CLINIC
#   }
  
#   required_cols <- c("sample", "pfs_code", "pfs_time", "disease_stage", "imwg")
#   if (!all(required_cols %in% colnames(clinic))) {
#     stop("Missing required clinical columns: sample, pfs_code, pfs_time, disease_stage, imwg")
#   }
  
#   ## Merge
#   all_data <- merge(clinic, matrix, by="sample")
  
#   # Merge CLASS if provided
#   if (!is.null(CLASS)) {
#     if (class(CLASS)=="character") {
#       class_df <- read.delim(CLASS, stringsAsFactors=FALSE)
#     } else {
#       class_df <- CLASS
#     }
#   }
  
#   if (!is.null(CLASS)) {
#     if (!"sample" %in% colnames(class_df)) stop("CLASS must contain a 'sample' column")
#     all_data <- merge(all_data, class_df, by = "sample", all.x = TRUE)
#   }
  
#   colnames(all_data)[ncol(all_data)]="class"
  
#   # Always apply tumor filtering first if CLASS is present
#   if (!is.null(CLASS) && "class" %in% colnames(all_data)) {
#     removed_classes <- table(all_data$class)
#     all_data <- all_data[all_data$class == "tumor", ]
#     cat("Tumor filtering applied.
# ")
#     cat("Samples removed by class filter:
# ")
#     print(removed_classes[names(removed_classes) != "tumor"])
#     cat("Remaining tumor samples:", nrow(all_data), "
# ")
#   }
#   #}  # Close tumor filtering and then CLASS
#   #}
#   #  ]
#   #}  # End of CLASS merge and tumor filtering
  
#   # Determine selection set
#   cat("
# === SELECTION FILTERING REPORT ===
# ")
#   initial_n <- nrow(all_data)
#   if (selection_by == "stage") {
#     data <- all_data[all_data$disease_stage %in% selection_value & !is.na(all_data$pfs_code), ]
#     cat("Selection by disease stage:", selection_value, "
# ")
#     cat("Samples removed by stage filter:", initial_n - nrow(data), "
# ")
#   } else if (selection_by == "cohort") {
#     data <- all_data[all_data$cohort %in% selection_value & !is.na(all_data$pfs_code), ]
#     cat("Selection by cohort:", selection_value, "
# ")
#     cat("Samples removed by cohort filter:", initial_n - nrow(data), "
# ")
#   } else {
#     data <- all_data[all_data$disease_stage %in% c("SMM", "MGUS") & !is.na(all_data$pfs_code), ]
#     cat("Default selection (SMM and MGUS).
# ")
#     cat("Samples removed by default stage filter:", initial_n - nrow(data), "
# ")
#   }
  
  
  
  
#   if (use_known_features) {
#     sig_features <- c("CNV.SNV_TET2", "CNV.SNV_NF1", "CNV_Del_10q24.32", "MYC", "CNV.Sig", "SNV_NRAS")
#     sig_features <- sig_features[sig_features %in% colnames(all_data)]
#     stat_df <- data.frame()
#     for (feature in sig_features) {
#       x <- data
#       x$feature <- x[[feature]]
#       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#         x$feature[x$feature > 1] <- 1
#       }
#       if (length(unique(x$feature)) < 2) next
#       fit <- coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       hr <- exp(coef(fit))
#       ci <- exp(confint(fit))
#       result_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       km_p_value <- as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
#       pval <- km_p_value
#       present <- sum(x$feature %in% c(1, 2), na.rm=TRUE)
#       absent <- sum(x$feature == 0, na.rm=TRUE)
#       present_MGUS <- sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
#       present_SMM <- sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
#       stat_df <- rbind(stat_df, data.frame(
#         feature = feature,
#         p_value = km_p_value,
#         HR = hr,
#         HR_low = ci[1],
#         HR_high = ci[2],
#         present = present,
#         absent = absent,
#         present_MGUS = present_MGUS,
#         present_SMM = present_SMM,
#         FDR = NA
#       ))
#     }
#     stat_df$FDR <- p.adjust(stat_df$p_value, method="fdr")
#     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
#   } else {
#     stat_df <- data.frame()
#     for (feature in colnames(data)[!colnames(data) %in% c("sample", "pfs_time", "pfs_code", "disease_stage", "cohort", "imwg", "class", "seq", "hy")]) {
#       x <- data
#       x$feature <- x[[feature]]
#       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#         x$feature[x$feature > 1] <- 1
#       }
#       if (length(unique(x$feature)) < 2) next
#       fit <- coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       hr <- exp(coef(fit))
#       ci <- exp(confint(fit))
#       result_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#       km_p_value <- as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
#       pval <- km_p_value
#       present <- sum(x$feature %in% c(1, 2), na.rm=TRUE)
#       absent <- sum(x$feature == 0, na.rm=TRUE)
#       present_MGUS <- sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
#       present_SMM <- sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
#       stat_df <- rbind(stat_df, data.frame(
#         feature = feature,
#         p_value = pval,
#         HR = hr,
#         HR_low = ci[1],
#         HR_high = ci[2],
#         present = present,
#         absent = absent,
#         present_MGUS = present_MGUS,
#         present_SMM = present_SMM,
#         FDR = NA
#       ))
#     }
#     stat_df$FDR <- p.adjust(stat_df$p_value, method = "fdr")
#     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
    
#     if (filter_by == "FDR") {
#       sig_features <- stat_df$feature[stat_df$FDR < filter_value & stat_df$present >= min_event_count]
#     } else {
#       sig_features <- stat_df$feature[stat_df$p_value < 0.05 & stat_df$present >= min_event_count]
#     }
#     sig_features <- sig_features[sig_features %in% colnames(all_data)]
#   }
  
#   if (length(sig_features) > 0) {
#     all_data$genomic_score <- rowSums(all_data[, sig_features, drop=FALSE], na.rm=TRUE)
#   } else {
#     warning("None of the selected features are present in the input data. Setting all genomic scores to 0.")
#     all_data$genomic_score <- 0
#   }
#   all_data$genomic_score[all_data$genomic_score > 0] <- 2
  
  
#   score_df <- data.frame(sample = all_data$sample,
#                          genomic_score = all_data$genomic_score,
#                          disease_stage = all_data$disease_stage,
#                          cohort = all_data$cohort,
#                          imwg = all_data$imwg)
#   write.table(score_df, file.path(output_dir, "genomic_score_per_sample.txt"), sep="	", row.names=FALSE, quote=FALSE)
  
#   # Barplot of HR
#   pdf(file.path(output_dir, "plots", "HR_barplot.pdf"), onefile = TRUE)
#   hr_data <- stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
#   if (nrow(hr_data) > 0) {
#     barplot(hr_data$HR,
#             names.arg = hr_data$feature,
#             las = 2, col = "coral", ylab = "Hazard Ratio")
#   } else {
#     plot.new()
#     text(0.5, 0.5, "No finite HR values with >=3 events to plot")
#   }
#   dev.off()
  
#   # Plot HR with CI error bars
#   pdf(file.path(output_dir, "plots", "HR_with_CI_plot.pdf"), onefile = TRUE)
#   stat_df_sig <- stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
#   if (nrow(stat_df_sig) > 0 && any(is.finite(stat_df_sig$HR_high))) {
#     par(mar = c(10, 5, 4, 2))
#     ylim_max <- max(stat_df_sig$HR_high, na.rm = TRUE)
#     if (is.finite(ylim_max)) {
#       x <- barplot(stat_df_sig$HR, names.arg = stat_df_sig$feature, las = 2,
#                    col = "skyblue", ylab = "Hazard Ratio", ylim = c(0, ylim_max * 1.1))
#       segments(x0 = x, x1 = x, y0 = stat_df_sig$HR_low, y1 = stat_df_sig$HR_high)
#     } else {
#       plot.new()
#       text(0.5, 0.5, "No finite HR CI to plot")
#     }
#   } else {
#     plot.new()
#     text(0.5, 0.5, "No significant features with valid HR and CI")
#   }
#   dev.off()
  
#   # KM plots for significant features
#   for (feature in sig_features) {
#     x <- data
#     x$feature <- x[[feature]]
#     if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
#       x$feature[x$feature > 1] <- 1
#     }
#     if (length(unique(x$feature)) < 2) next
#     fit_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
#     pdf(file.path(output_dir, "km_plots", paste0("KM_", gsub("/", "_", feature), ".pdf")), onefile = FALSE)
#     print(ggsurvplot(fit_km, data = x, pval = TRUE, risk.table = TRUE, xlab = "Time", ylab = feature))
#     dev.off()
    
    
#   }
  
#   # Multivariate Cox model using all samples
#   if (length(sig_features) >= 1 && any(!is.na(all_data$pfs_code))) {
#     full_data <- all_data[!is.na(all_data$pfs_code), ]
#     if (nrow(full_data) > 0) {
#       cox_all <- coxph(Surv(pfs_time, pfs_code) ~ genomic_score + imwg, data = full_data, robust = TRUE)
#       summary_all <- as.data.frame(summary(cox_all)$coefficients)
#       summary_all$variable <- rownames(summary_all)
#       write.table(summary_all, file.path(output_dir, "cox_analysis", "cox_summary_all_samples.txt"), sep="	", row.names=FALSE, quote=FALSE)
      
#       pdf(file.path(output_dir, "cox_analysis", "forest_plot_all_samples.pdf"))
#       print(forest_model(cox_all))
#       dev.off()
#     }
#   }
  
#   result <- list(
#     score_df = score_df,
#     analysis_data = data
#   )
#   return(result)
# }



# # genomic_score <- function(MATRIX, CLINIC, CLASS = NULL, output_dir = "genomic_score_output", filter_value = 0.1,
# #                           selection_by = "all", selection_value = NULL,
# #                           filter_by = "FDR", min_event_count = 2, use_known_features = FALSE) {
# #   # Required packages
# #   required_packages <- c("survival", "forestmodel", "survminer")
# #   for (pkg in required_packages) {
# #     if (!requireNamespace(pkg, quietly = TRUE)) {
# #       install.packages(pkg)
# #     }
# #     library(pkg, character.only = TRUE)
# #   }
# #   
# #   # Create output subdirectories if they don't exist
# #   dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
# #   dir.create(file.path(output_dir, "km_plots"), showWarnings = FALSE, recursive = TRUE)
# #   dir.create(file.path(output_dir, "cox_analysis"), showWarnings = FALSE, recursive = TRUE)
# #   
# #   if (class(MATRIX)=="character") {
# #     matrix <- read.delim(MATRIX, stringsAsFactors=FALSE)
# #   } else {
# #     matrix <- MATRIX
# #   }
# #   
# #   if (class(CLINIC)=="character") {
# #     clinic <- read.delim(CLINIC, stringsAsFactors=FALSE)
# #   } else {
# #     clinic <- CLINIC
# #   }
# #   
# #   required_cols <- c("sample", "pfs_code", "pfs_time", "disease_stage", "imwg")
# #   if (!all(required_cols %in% colnames(clinic))) {
# #     stop("Missing required clinical columns: sample, pfs_code, pfs_time, disease_stage, imwg")
# #   }
# #   
# #   ## Merge
# #   all_data <- merge(clinic, matrix, by="sample")
# #   
# #   # Merge CLASS if provided
# #   if (!is.null(CLASS)) {
# #     if (class(CLASS)=="character") {
# #       class_df <- read.delim(CLASS, stringsAsFactors=FALSE)
# #     } else {
# #       class_df <- CLASS
# #     }
# #   }
# #   
# #   if (!is.null(CLASS)) {
# #     if (!"sample" %in% colnames(class_df)) stop("CLASS must contain a 'sample' column")
# #     all_data <- merge(all_data, class_df, by = "sample", all.x = TRUE)
# #   }
# #   
# #   colnames(all_data)[ncol(all_data)]="class"
# #   
# #   # Always apply tumor filtering first if CLASS is present
# #   if (!is.null(CLASS) && "class" %in% colnames(all_data)) {
# #     removed_classes <- table(all_data$class)
# #     all_data <- all_data[all_data$class == "tumor", ]
# #     cat("Tumor filtering applied.
# # ")
# #     cat("Samples removed by class filter:
# # ")
# #     print(removed_classes[names(removed_classes) != "tumor"])
# #     cat("Remaining tumor samples:", nrow(all_data), "
# # ")
# #   }
# #   #}  # Close tumor filtering and then CLASS
# #   #}
# #   #  ]
# #   #}  # End of CLASS merge and tumor filtering
# #   
# #   # Determine selection set
# #   cat("
# # === SELECTION FILTERING REPORT ===
# # ")
# #   initial_n <- nrow(all_data)
# #   if (selection_by == "stage") {
# #     data <- all_data[all_data$disease_stage %in% selection_value & !is.na(all_data$pfs_code), ]
# #     cat("Selection by disease stage:", selection_value, "
# # ")
# #     cat("Samples removed by stage filter:", initial_n - nrow(data), "
# # ")
# #   } else if (selection_by == "cohort") {
# #     data <- all_data[all_data$cohort %in% selection_value & !is.na(all_data$pfs_code), ]
# #     cat("Selection by cohort:", selection_value, "
# # ")
# #     cat("Samples removed by cohort filter:", initial_n - nrow(data), "
# # ")
# #   } else {
# #     data <- all_data[all_data$disease_stage %in% c("SMM", "MGUS") & !is.na(all_data$pfs_code), ]
# #     cat("Default selection (SMM and MGUS).
# # ")
# #     cat("Samples removed by default stage filter:", initial_n - nrow(data), "
# # ")
# #   }
# #   
# #   
# #   
# #   
# #   if (use_known_features) {
# #     sig_features <- c("CNV.SNV_TET2", "CNV.SNV_NF1", "CNV_Del_10q24.32", "MYC", "CNV.Sig", "SNV_NRAS")
# #     sig_features <- sig_features[sig_features %in% colnames(all_data)]
# #     stat_df <- data.frame()
# #     for (feature in sig_features) {
# #       x <- data
# #       x$feature <- x[[feature]]
# #       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
# #         x$feature[x$feature > 1] <- 1
# #       }
# #       if (length(unique(x$feature)) < 2) next
# #       fit <- coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
# #       hr <- exp(coef(fit))
# #       ci <- exp(confint(fit))
# #       result_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
# #       km_p_value <- as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
# #       pval <- km_p_value
# #       present <- sum(x$feature %in% c(1, 2), na.rm=TRUE)
# #       absent <- sum(x$feature == 0, na.rm=TRUE)
# #       present_MGUS <- sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
# #       present_SMM <- sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
# #       stat_df <- rbind(stat_df, data.frame(
# #         feature = feature,
# #         p_value = km_p_value,
# #         HR = hr,
# #         HR_low = ci[1],
# #         HR_high = ci[2],
# #         present = present,
# #         absent = absent,
# #         present_MGUS = present_MGUS,
# #         present_SMM = present_SMM,
# #         FDR = NA
# #       ))
# #     }
# #     stat_df$FDR <- p.adjust(stat_df$p_value, method="fdr")
# #     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
# #   } else {
# #     stat_df <- data.frame()
# #     for (feature in colnames(data)[!colnames(data) %in% c("sample", "pfs_time", "pfs_code", "disease_stage", "cohort", "imwg", "class", "seq", "hy")]) {
# #       x <- data
# #       x$feature <- x[[feature]]
# #       if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
# #         x$feature[x$feature > 1] <- 1
# #       }
# #       if (length(unique(x$feature)) < 2) next
# #       fit <- coxph(Surv(pfs_time, pfs_code) ~ feature, data = x)
# #       hr <- exp(coef(fit))
# #       ci <- exp(confint(fit))
# #       result_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
# #       km_p_value <- as.numeric(gsub("p = ", "", surv_pvalue(result_km)[,2]))
# #       pval <- km_p_value
# #       present <- sum(x$feature %in% c(1, 2), na.rm=TRUE)
# #       absent <- sum(x$feature == 0, na.rm=TRUE)
# #       present_MGUS <- sum(x$feature %in% c(1,2) & x$disease_stage == "MGUS", na.rm=TRUE)
# #       present_SMM <- sum(x$feature %in% c(1,2) & x$disease_stage == "SMM", na.rm=TRUE)
# #       stat_df <- rbind(stat_df, data.frame(
# #         feature = feature,
# #         p_value = pval,
# #         HR = hr,
# #         HR_low = ci[1],
# #         HR_high = ci[2],
# #         present = present,
# #         absent = absent,
# #         present_MGUS = present_MGUS,
# #         present_SMM = present_SMM,
# #         FDR = NA
# #       ))
# #     }
# #     stat_df$FDR <- p.adjust(stat_df$p_value, method = "fdr")
# #     write.table(stat_df, file.path(output_dir, "genomic_feature_stats.txt"), sep="	", row.names=FALSE, quote=FALSE)
# #     
# #     if (filter_by == "FDR") {
# #       sig_features <- stat_df$feature[stat_df$FDR < filter_value & stat_df$present >= min_event_count]
# #     } else {
# #       sig_features <- stat_df$feature[stat_df$p_value < 0.05 & stat_df$present >= min_event_count]
# #     }
# #     sig_features <- sig_features[sig_features %in% colnames(all_data)]
# #   }
# #   
# #   if (length(sig_features) > 0) {
# #     all_data$genomic_score <- rowSums(all_data[, sig_features, drop=FALSE], na.rm=TRUE)
# #   } else {
# #     warning("None of the selected features are present in the input data. Setting all genomic scores to 0.")
# #     all_data$genomic_score <- 0
# #   }
# #   all_data$genomic_score[all_data$genomic_score > 0] <- 2
# #   
# #   
# #   score_df <- data.frame(sample = all_data$sample,
# #                          genomic_score = all_data$genomic_score,
# #                          disease_stage = all_data$disease_stage,
# #                          cohort = all_data$cohort,
# #                          imwg = all_data$imwg)
# #   write.table(score_df, file.path(output_dir, "genomic_score_per_sample.txt"), sep="	", row.names=FALSE, quote=FALSE)
# #   
# #   # Barplot of HR
# #   pdf(file.path(output_dir, "plots", "HR_barplot.pdf"), onefile = TRUE)
# #   hr_data <- stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
# #   if (nrow(hr_data) > 0) {
# #     barplot(hr_data$HR,
# #             names.arg = hr_data$feature,
# #             las = 2, col = "coral", ylab = "Hazard Ratio")
# #   } else {
# #     plot.new()
# #     text(0.5, 0.5, "No finite HR values with >=3 events to plot")
# #   }
# #   dev.off()
# #   
# #   # Plot HR with CI error bars
# #   pdf(file.path(output_dir, "plots", "HR_with_CI_plot.pdf"), onefile = TRUE)
# #   stat_df_sig <- stat_df[stat_df$FDR < filter_value & stat_df$present >= 3 & is.finite(stat_df$HR), ]
# #   if (nrow(stat_df_sig) > 0 && any(is.finite(stat_df_sig$HR_high))) {
# #     par(mar = c(10, 5, 4, 2))
# #     ylim_max <- max(stat_df_sig$HR_high, na.rm = TRUE)
# #     if (is.finite(ylim_max)) {
# #       x <- barplot(stat_df_sig$HR, names.arg = stat_df_sig$feature, las = 2,
# #                    col = "skyblue", ylab = "Hazard Ratio", ylim = c(0, ylim_max * 1.1))
# #       segments(x0 = x, x1 = x, y0 = stat_df_sig$HR_low, y1 = stat_df_sig$HR_high)
# #     } else {
# #       plot.new()
# #       text(0.5, 0.5, "No finite HR CI to plot")
# #     }
# #   } else {
# #     plot.new()
# #     text(0.5, 0.5, "No significant features with valid HR and CI")
# #   }
# #   dev.off()
# #   
# #   # KM plots for significant features
# #   for (feature in sig_features) {
# #     x <- data
# #     x$feature <- x[[feature]]
# #     if (length(unique(x$feature)) == 3 && all(unique(x$feature) %in% c(0,1,2))) {
# #       x$feature[x$feature > 1] <- 1
# #     }
# #     if (length(unique(x$feature)) < 2) next
# #     fit_km <- surv_fit(Surv(pfs_time, pfs_code) ~ feature, data = x)
# #     pdf(file.path(output_dir, "km_plots", paste0("KM_", gsub("/", "_", feature), ".pdf")), onefile = FALSE)
# #     print(ggsurvplot(fit_km, data = x, pval = TRUE, risk.table = TRUE, xlab = "Time", ylab = feature))
# #     dev.off()
# #     
# #     
# #   }
# #   
# #   # Multivariate Cox model using all samples
# #   if (length(sig_features) >= 1 && any(!is.na(all_data$pfs_code))) {
# #     full_data <- all_data[!is.na(all_data$pfs_code), ]
# #     if (nrow(full_data) > 0) {
# #       cox_all <- coxph(Surv(pfs_time, pfs_code) ~ genomic_score + imwg, data = full_data, robust = TRUE)
# #       summary_all <- as.data.frame(summary(cox_all)$coefficients)
# #       summary_all$variable <- rownames(summary_all)
# #       write.table(summary_all, file.path(output_dir, "cox_analysis", "cox_summary_all_samples.txt"), sep="	", row.names=FALSE, quote=FALSE)
# #       
# #       pdf(file.path(output_dir, "cox_analysis", "forest_plot_all_samples.pdf"))
# #       print(forest_model(cox_all))
# #       dev.off()
# #     }
# #   }
# #   
# #   return(score_df)
# # }
