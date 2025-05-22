###########################################################################
############################################################
#############################################
##############################
###############
#### . . . . - Tumor/No.Tumor classification - . . . . ####
###############
##############################
#############################################
############################################################
###########################################################################

# . . . . . . . . . . -> -> -> -> -> add the method to compute the significant features

SMM_TnT_Cla=function(GC_input_matrix,cyto,sig.features=NULL,assembly=NULL) {
  
  output=list()
  
  #### . . . - Preparation Data - . . . ####
  if(class(GC_input_matrix)=="character"){
    matrix=read.delim(GC_input_matrix,
                      stringsAsFactors= F, header=T)
  }else{
    matrix=GC_input_matrix
  }
  
  if(ncol(matrix) != 133){
    stop("ERROR: Please provide the input matrix with all 138 data columns")
  }
  
  if(class(cyto)=="character"){
    new.features=read.delim(cyto,
                            stringsAsFactors= F, header=T)
  }else{
    new.features=cyto
  }
  
  if(sum(c("sample", "CCND3", "CCND2", "MAFA", "MYC","hy") %in% names(new.features)) != 6){
    stop("ERROR: Please provide the following data columns: sample, CCND3, CCND2, MAFA, MYC, hy")
  }
  
  if(!is.null(sig.features)){
    sig_feat=sig.features
  }else if(assembly %in% c("hg38","HG38","Hg38","GRCh38","GRCH38")){
    # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
    #            "CNV_Gain_Amp1q",
    #            "CNV_chr8p.loss",
    #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
    #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
    #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
    #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_H1.3",
    #            "APOBEC",
    #            "t_MMSET",
    #            "CNV.Sig",
    #            "hy")
    # 2       sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
    #            "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
    #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
    #            "SNV_FGFR3","SNV_NRAS",
    #            "CNV.Sig",
    #            "APOBEC",
    #            "t_MMSET","MYC",
    #            "hy")
      sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_6q26","CNV_Del_8q24.21", 
               "CNV.SNV_ARID2","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_TET2","SNV_FGFR3",       
               "SNV_NRAS","CNV.Sig","APOBEC","t_MMSET","MYC","hy")
  }else{
    # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
    #            "CNV_Gain_Amp1q",
    #            "CNV_chr8p.loss",
    #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
    #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
    #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
    #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
    #            "APOBEC",
    #            "t_MMSET",
    #            "CNV.Sig",
    #            "hy")
    # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
    #            "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
    #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
    #            "SNV_FGFR3","SNV_NRAS",
    #            "CNV.Sig",
    #            "APOBEC",
    #            "t_MMSET","MYC",
    #            "hy")
    sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain","CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_6q26","CNV_Del_8q24.21", 
      "CNV.SNV_ARID2","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_TET2","SNV_FGFR3",       
      "SNV_NRAS","CNV.Sig","APOBEC","t_MMSET","MYC","hy")
  }
  
  # if(assembly %in% c("hg38","HG38","Hg38","GRCh38","GRCH38")){
  #   # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
  #   #            "CNV_Gain_Amp1q",
  #   #            "CNV_chr8p.loss",
  #   #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
  #   #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
  #   #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
  #   #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_H1.3",
  #   #            "APOBEC",
  #   #            "t_MMSET",
  #   #            "CNV.Sig",
  #   #            "hy")
  #   sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
  #              "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
  #              "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
  #              "SNV_FGFR3","SNV_NRAS",
  #              "CNV.Sig",
  #              "APOBEC",
  #              "t_MMSET","MYC",
  #              "hy")
  #   
  #   
  # }else{
  #   sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
  #              "CNV_Gain_Amp1q",
  #              "CNV_chr8p.loss",
  #              "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
  #              "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
  #              "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
  #              "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
  #              "APOBEC",
  #              "t_MMSET",
  #              "CNV.Sig",
  #              "hy")
  # }
  
  #### . . . - Merging all the information - . . . ####
  matrix=merge(matrix,new.features,by="sample",all.x=T)
  dim(matrix)
  
  #### . . . - Creating the genomicd object - . . . ####
  data=matrix
  copynumbers=colnames(data)[grep("CNV_",colnames(data))]
  tumorsuppressorgenes=colnames(data)[grep("CNV.SNV_",colnames(data))]
  oncogenesothers=c(colnames(data)[grep("^SNV_|^SNV.",colnames(data))])
  genomics=data[,c(colnames(data)[grep("sample|sampleID|Sample|SampleID|SAMPLE|SAMPLE.ID",colnames(data))],
                   copynumbers,tumorsuppressorgenes,oncogenesothers,
                   "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy")]
  
  
  
  ######################################
  #### . . . - Starting the tumor/no.tumor classification - . . . ####
  ######################################
  #### . . . - Defining tumor: considering the presence of: MYC, APOBEC, CNV.Sig, KRAS/NRAS/BRAF mutations, t_MMSET and TP53 - . . . ####
  genomics$class.new=NA
  genomics$class.new[genomics$MYC %in% c(1)]="tumor"
  genomics$class.new[genomics$APOBEC %in% c(1)]="tumor"
  genomics$class.new[genomics$CNV.Sig %in% c(1)]="tumor" 
  genomics$class.new[(genomics$SNV_KRAS +
                        genomics$SNV_NRAS +
                        genomics$SNV_BRAF) > 0]="tumor"
  genomics$class.new[genomics$t_MMSET %in% c(1)]="tumor"
  genomics$class.new[genomics$CNV.SNV_TP53 > 0]="tumor"
  table(genomics$class.new,useNA="ifany")
  # tumor    NA
  # 253     139
  
  #### . . . - Initializing a new object with the sample and class information - . . . ####
  summary.file.class=genomics[,c("sample","class.new")] ## reference file
  
  #### . . . - Defining tumor by counting the DRIVERs - . . . ####
  genomics_check=genomics[is.na(genomics$class.new),]
  
  for(i in 2:138){
    genomics_check[,i]=as.numeric(as.character(genomics_check[,i]))
    genomics_check[,i]=ifelse(is.na(genomics_check[,i]),0,genomics_check[,i])
  }
  
  #### . . . - Selecting the non-canonical drivers (without sampleID column) . . . - ####
  non.canonical.drivers=setdiff(colnames(genomics_check)[-c(1,139)],sig_feat)
  
  # #### . . . - Creating all the possible variables to finish the classification . . . - ####
  # genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-37]]) # exclude HY
  # genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr21.gain")])# include HY and exclude all gains
  # genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 37)]])# exclude HY and exclude all gains except for +18
  # genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
  # genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
  # genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
  # genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
  #   genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  
  # #### . . . - Creating all the possible variables to finish the classification . . . - ####
  # genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-27]]) # exclude HY
  # genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr11.gain")])# include HY and exclude all gains
  # genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 27)]])# exclude HY and exclude all gains except for +18
  # genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
  # genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
  # genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
  # genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
  #   genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  # #### . . . - Creating all the possible variables to finish the classification . . . - ####
  # genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-28]]) # exclude HY
  # genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr7.gain")])# include HY and exclude all gains
  # genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 28)]])# exclude HY and exclude all gains except for +18
  # genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
  # genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
  # genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
  # genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
  #   genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  #### May 22 2025
  #### . . . - Creating all the possible variables to finish the classification . . . - ####
  genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-28]]) # exclude HY
  genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr7.gain")])# include HY and exclude all gains
  genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 28)]])# exclude HY ####and exclude all gains except for +18
  genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,105:109)])]) # chr 18 included as well as all non HY gains
  genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(10:104)])]) ###chr##18##included###
  genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
  genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
    genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  
  #### . . . - Defining tumor: considering the presence of: if you have feature associated with SMM progression + hy or + igh translocation you are a tumor - . . . ####
  #dim(genomics_check) # 139
  genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$igh > 0]="tumor" # IGH + at.least.1.driver
  genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$hy > 0]="tumor"  # HY + at.least.1.driver 
  genomics_check$class.new[genomics_check$hy > 0 & genomics_check$igh > 0]="tumor"                      # HY + IGH
  
  #### . . . - Defining tumor: considering the presence of: if you have at least one MM prognostic features + any other feature you are tumor  - . . . ####
  # . - 1 is because the prognostic driver is counted twice (BAC: correct)
  genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$non_canonical_driver > 0]="tumor"
  #table(summary.file.class$class.new,useNA="ifany")
  summary.file.class$class.new[summary.file.class$sample %in% genomics_check$sample[genomics_check$class.new=="tumor"]]="tumor"
  #table(summary.file.class$class.new)
  # 297 tumor
  
  #### . . . - Defining non_tumor: considering the presence of: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
  # . - Define as non tumor all the samples with no drivers, but passed by filters
  genomics_check.2=genomics_check[is.na(genomics_check$class.new),]
  #dim(genomics_check.2)
  genomics_check.2$class.new[genomics_check.2$driver_count_all_gains == 0 &
                               genomics_check.2$driver_count_only_hy == 0 &
                               genomics_check.2$driver_count_only_no_hy == 0 &
                               genomics_check.2$non_canonical_driver == 0 &
                               genomics_check.2$non_canonical_driver_no_gains &
                               genomics_check.2$igh == 0]="non_tumor"
  
  
  summary.file.class$class.new[summary.file.class$sample %in% genomics_check.2$class.new]="non_tumor"
  # table(summary.file.class$class.new,useNA="ifany")
  # tumor  <NA> 
  #   297    95
  
  #### . . . - Defining pts: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
  genomics_check.3=genomics_check.2[is.na(genomics_check.2$class),]
  #dim(genomics_check.3) #95
  # genomics_check.3[genomics_check.3$driver_count_only_no_hy > 0,] ## all the samples with at least one genomic driver associated with progression are tumor
  genomics_check.3$non_canonical_driver_no_gains2=genomics_check.3$non_canonical_driver_no_gains - genomics_check.3$CNV.SNV_TGDS # we consider only one event in 13q
  #genomics_check.3[genomics_check.3$igh > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample # or >2
  #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample #or >2
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
                                                                                               genomics_check.3$non_canonical_driver > 1))]="tumor"
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
                                                                                                genomics_check.3$non_canonical_driver > 1))]="tumor"
  #genomics_check.3[genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0,]$sample
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0)]="tumor"
  
  ##### presence of >2 non canonical is tumor - sse above
  #### presence of 2 + HY and IGH 
  
  #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2),]$sample
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2))]="tumor" 
  #table(genomics_check.3$class.new,useNA="ifany")
  # tumor  <NA> 
  #   15    80
  
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains %in% c(0,1,2))]="non_tumor"  #c(0,1,2)
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(0,1,2))]="non_tumor" #c(0,1,2)
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
  genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
  #table(genomics_check.3$class.new,useNA="ifany")
  # non_tumor     tumor     
  #        80        15      
  #genomics_check.3[is.na(genomics_check.3$class.new),]
  
  summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="tumor" ]]="tumor"
  summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="non_tumor" ]]="non_tumor"
  
  table(summary.file.class$class.new,useNA="ifany")
  # non_tumor     tumor 
  #        80       312 
  
  output$class=summary.file.class
  output$matrix=genomics[,-ncol(genomics)]
  
  return(output)
}























# # . . . . . . . . . . -> -> -> -> -> add the method to compute the significant features

# SMM_TnT_Cla=function(GC_input_matrix,cyto,sig.features=NULL,assembly=NULL) {
  
#   output=list()
  
#   #### . . . - Preparation Data - . . . ####
#   if(class(GC_input_matrix)=="character"){
#     matrix=read.delim(GC_input_matrix,
#                       stringsAsFactors= F, header=T)
#   }else{
#     matrix=GC_input_matrix
#   }
  
#   if(ncol(matrix) != 133){
#     stop("ERROR: Please provide the input matrix with all 138 data columns")
#   }
  
#   if(class(cyto)=="character"){
#     new.features=read.delim(cyto,
#                             stringsAsFactors= F, header=T)
#   }else{
#     new.features=cyto
#   }
  
#   if(sum(c("sample", "CCND3", "CCND2", "MAFA", "MYC","hy") %in% names(new.features)) != 6){
#     stop("ERROR: Please provide the following data columns: sample, CCND3, CCND2, MAFA, MYC, hy")
#   }
  
#   if(!is.null(sig.features)){
#     sig_feat=sig.features
#   }else if(assembly %in% c("hg38","HG38","Hg38","GRCh38","GRCH38")){
#     # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
#     #            "CNV_Gain_Amp1q",
#     #            "CNV_chr8p.loss",
#     #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
#     #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
#     #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
#     #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_H1.3",
#     #            "APOBEC",
#     #            "t_MMSET",
#     #            "CNV.Sig",
#     #            "hy")
#     sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
#                "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
#                "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
#                "SNV_FGFR3","SNV_NRAS",
#                "CNV.Sig",
#                "APOBEC",
#                "t_MMSET","MYC",
#                "hy")
#   }else{
#     # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
#     #            "CNV_Gain_Amp1q",
#     #            "CNV_chr8p.loss",
#     #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
#     #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
#     #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
#     #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
#     #            "APOBEC",
#     #            "t_MMSET",
#     #            "CNV.Sig",
#     #            "hy")
#     sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
#                "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
#                "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
#                "SNV_FGFR3","SNV_NRAS",
#                "CNV.Sig",
#                "APOBEC",
#                "t_MMSET","MYC",
#                "hy")
#   }
  
#   # if(assembly %in% c("hg38","HG38","Hg38","GRCh38","GRCH38")){
#   #   # sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
#   #   #            "CNV_Gain_Amp1q",
#   #   #            "CNV_chr8p.loss",
#   #   #            "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
#   #   #            "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
#   #   #            "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
#   #   #            "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_H1.3",
#   #   #            "APOBEC",
#   #   #            "t_MMSET",
#   #   #            "CNV.Sig",
#   #   #            "hy")
#   #   sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr21.gain",  
#   #              "CNV_Del_10q24.32","CNV_Del_2q31.1","CNV_Del_8q24.21",
#   #              "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_CREBBP","CNV.SNV_CYLD","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_TET2",
#   #              "SNV_FGFR3","SNV_NRAS",
#   #              "CNV.Sig",
#   #              "APOBEC",
#   #              "t_MMSET","MYC",
#   #              "hy")
#   #   
#   #   
#   # }else{
#   #   sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
#   #              "CNV_Gain_Amp1q",
#   #              "CNV_chr8p.loss",
#   #              "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
#   #              "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
#   #              "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
#   #              "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
#   #              "APOBEC",
#   #              "t_MMSET",
#   #              "CNV.Sig",
#   #              "hy")
#   # }
  
#   #### . . . - Merging all the information - . . . ####
#   matrix=merge(matrix,new.features,by="sample",all.x=T)
#   dim(matrix)
  
#   #### . . . - Creating the genomicd object - . . . ####
#   data=matrix
#   copynumbers=colnames(data)[grep("CNV_",colnames(data))]
#   tumorsuppressorgenes=colnames(data)[grep("CNV.SNV_",colnames(data))]
#   oncogenesothers=c(colnames(data)[grep("^SNV_|^SNV.",colnames(data))])
#   genomics=data[,c(colnames(data)[grep("sample|sampleID|Sample|SampleID|SAMPLE|SAMPLE.ID",colnames(data))],
#                    copynumbers,tumorsuppressorgenes,oncogenesothers,
#                    "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy")]
  
  
  
#   ######################################
#   #### . . . - Starting the tumor/no.tumor classification - . . . ####
#   ######################################
#   #### . . . - Defining tumor: considering the presence of: MYC, APOBEC, CNV.Sig, KRAS/NRAS/BRAF mutations, t_MMSET and TP53 - . . . ####
#   genomics$class.new=NA
#   genomics$class.new[genomics$MYC %in% c(1)]="tumor"
#   genomics$class.new[genomics$APOBEC %in% c(1)]="tumor"
#   genomics$class.new[genomics$CNV.Sig %in% c(1)]="tumor" 
#   genomics$class.new[(genomics$SNV_KRAS +
#                         genomics$SNV_NRAS +
#                         genomics$SNV_BRAF) > 0]="tumor"
#   genomics$class.new[genomics$t_MMSET %in% c(1)]="tumor"
#   genomics$class.new[genomics$CNV.SNV_TP53 > 0]="tumor"
#   table(genomics$class.new,useNA="ifany")
#   # tumor    NA
#   # 253     139
  
#   #### . . . - Initializing a new object with the sample and class information - . . . ####
#   summary.file.class=genomics[,c("sample","class.new")] ## reference file
  
#   #### . . . - Defining tumor by counting the DRIVERs - . . . ####
#   genomics_check=genomics[is.na(genomics$class.new),]
  
#   for(i in 2:138){
#     genomics_check[,i]=as.numeric(as.character(genomics_check[,i]))
#     genomics_check[,i]=ifelse(is.na(genomics_check[,i]),0,genomics_check[,i])
#   }
  
#   #### . . . - Selecting the non-canonical drivers (without sampleID column) . . . - ####
#   non.canonical.drivers=setdiff(colnames(genomics_check)[-c(1,139)],sig_feat)
  
#   # #### . . . - Creating all the possible variables to finish the classification . . . - ####
#   # genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-37]]) # exclude HY
#   # genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr21.gain")])# include HY and exclude all gains
#   # genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 37)]])# exclude HY and exclude all gains except for +18
#   # genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
#   # genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
#   # genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
#   # genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
#   #   genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  
#   #### . . . - Creating all the possible variables to finish the classification . . . - ####
#   genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-27]]) # exclude HY
#   genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr11.gain")])# include HY and exclude all gains
#   genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 27)]])# exclude HY and exclude all gains except for +18
#   genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
#   genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
#   genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
#   genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
#     genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
  
#   #### . . . - Defining tumor: considering the presence of: if you have feature associated with SMM progression + hy or + igh translocation you are a tumor - . . . ####
#   #dim(genomics_check) # 139
#   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$igh > 0]="tumor" # IGH + at.least.1.driver
#   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$hy > 0]="tumor"  # HY + at.least.1.driver 
#   genomics_check$class.new[genomics_check$hy > 0 & genomics_check$igh > 0]="tumor"                      # HY + IGH
  
#   #### . . . - Defining tumor: considering the presence of: if you have at least one MM prognostic features + any other feature you are tumor  - . . . ####
#   # . - 1 is because the prognostic driver is counted twice (BAC: correct)
#   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$non_canonical_driver > 0]="tumor"
#   #table(summary.file.class$class.new,useNA="ifany")
#   summary.file.class$class.new[summary.file.class$sample %in% genomics_check$sample[genomics_check$class.new=="tumor"]]="tumor"
#   #table(summary.file.class$class.new)
#   # 297 tumor
  
#   #### . . . - Defining non_tumor: considering the presence of: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
#   # . - Define as non tumor all the samples with no drivers, but passed by filters
#   genomics_check.2=genomics_check[is.na(genomics_check$class.new),]
#   #dim(genomics_check.2)
#   genomics_check.2$class.new[genomics_check.2$driver_count_all_gains == 0 &
#                                genomics_check.2$driver_count_only_hy == 0 &
#                                genomics_check.2$driver_count_only_no_hy == 0 &
#                                genomics_check.2$non_canonical_driver == 0 &
#                                genomics_check.2$non_canonical_driver_no_gains &
#                                genomics_check.2$igh == 0]="non_tumor"
  
  
#   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.2$class.new]="non_tumor"
#   # table(summary.file.class$class.new,useNA="ifany")
#   # tumor  <NA> 
#   #   297    95
  
#   #### . . . - Defining pts: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
#   genomics_check.3=genomics_check.2[is.na(genomics_check.2$class),]
#   #dim(genomics_check.3) #95
#   # genomics_check.3[genomics_check.3$driver_count_only_no_hy > 0,] ## all the samples with at least one genomic driver associated with progression are tumor
#   genomics_check.3$non_canonical_driver_no_gains2=genomics_check.3$non_canonical_driver_no_gains - genomics_check.3$CNV.SNV_TGDS # we consider only one event in 13q
#   #genomics_check.3[genomics_check.3$igh > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample # or >2
#   #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample #or >2
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
#                                                                                                genomics_check.3$non_canonical_driver > 1))]="tumor"
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
#                                                                                                 genomics_check.3$non_canonical_driver > 1))]="tumor"
#   #genomics_check.3[genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0,]$sample
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$non_canonical_driver_no_gains2 > 1 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0)]="tumor"
  
#   ##### presence of >2 non canonical is tumor - sse above
#   #### presence of 2 + HY and IGH 
  
#   #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2),]$sample
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2))]="tumor" 
#   #table(genomics_check.3$class.new,useNA="ifany")
#   # tumor  <NA> 
#   #   15    80
  
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains %in% c(0,1,2))]="non_tumor"  #c(0,1,2)
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(0,1,2))]="non_tumor" #c(0,1,2)
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
#   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
#   #table(genomics_check.3$class.new,useNA="ifany")
#   # non_tumor     tumor     
#   #        80        15      
#   #genomics_check.3[is.na(genomics_check.3$class.new),]
  
#   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="tumor" ]]="tumor"
#   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="non_tumor" ]]="non_tumor"
  
#   table(summary.file.class$class.new,useNA="ifany")
#   # non_tumor     tumor 
#   #        80       312 
  
#   output$class=summary.file.class
#   output$matrix=genomics[,-ncol(genomics)]
  
#   return(output)
# }
















# # SMM_TnT_Cla=function(GC_input_matrix,cyto,assembly) {
  
# #   output=list()
  
# #   #### . . . - Preparation Data - . . . ####
# #   if(class(GC_input_matrix)=="character"){
# #     matrix=read.delim(GC_input_matrix,
# #                       stringsAsFactors= F, header=T)
# #   }else{
# #     matrix=GC_input_matrix
# #   }
  
# #   if(ncol(matrix) != 133){
# #     stop("ERROR: Please provide the input matrix with all 138 data columns")
# #   }
  
# #   if(class(cyto)=="character"){
# #     new.features=read.delim(cyto,
# #                             stringsAsFactors= F, header=T)
# #   }else{
# #     new.features=cyto
# #   }
  
# #   if(sum(c("sample", "CCND3", "CCND2", "MAFA", "MYC","hy") %in% names(new.features)) != 6){
# #     stop("ERROR: Please provide the following data columns: sample, CCND3, CCND2, MAFA, MYC, hy")
# #   }
  
# #   if(assembly %in% c("hg38","HG38","Hg38","GRCh38","GRCH38")){
# #     sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
# #                "CNV_Gain_Amp1q",
# #                "CNV_chr8p.loss",
# #                "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
# #                "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
# #                "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
# #                "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_H1.3",
# #                "APOBEC",
# #                "t_MMSET",
# #                "CNV.Sig",
# #                "hy")
# #   }else{
# #     sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
# #                "CNV_Gain_Amp1q",
# #                "CNV_chr8p.loss",
# #                "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
# #                "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_FAM46C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1",
# #                "CNV.SNV_TET2","CNV.SNV_TP53","CNV.SNV_ZNF292",
# #                "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
# #                "APOBEC",
# #                "t_MMSET",
# #                "CNV.Sig",
# #                "hy")
# #   }
  
# #   #### . . . - Merging all the information - . . . ####
# #   matrix=merge(matrix,new.features,by="sample",all.x=T)
# #   dim(matrix)
  
# #   #### . . . - Creating the genomicd object - . . . ####
# #   data=matrix
# #   copynumbers=colnames(data)[grep("CNV_",colnames(data))]
# #   tumorsuppressorgenes=colnames(data)[grep("CNV.SNV_",colnames(data))]
# #   oncogenesothers=c(colnames(data)[grep("^SNV_|^SNV.",colnames(data))])
# #   genomics=data[,c(colnames(data)[grep("sample|sampleID|Sample|SampleID|SAMPLE|SAMPLE.ID",colnames(data))],
# #                    copynumbers,tumorsuppressorgenes,oncogenesothers,
# #                    "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy")]
  
  
  
# #   ######################################
# #   #### . . . - Starting the tumor/no.tumor classification - . . . ####
# #   ######################################
# #   #### . . . - Defining tumor: considering the presence of: MYC, APOBEC, CNV.Sig, KRAS/NRAS/BRAF mutations, t_MMSET and TP53 - . . . ####
# #   genomics$class.new=NA
# #   genomics$class.new[genomics$MYC %in% c(1)]="tumor"
# #   genomics$class.new[genomics$APOBEC %in% c(1)]="tumor"
# #   genomics$class.new[genomics$CNV.Sig %in% c(1)]="tumor" 
# #   genomics$class.new[(genomics$SNV_KRAS +
# #                         genomics$SNV_NRAS +
# #                         genomics$SNV_BRAF) > 0]="tumor"
# #   genomics$class.new[genomics$t_MMSET %in% c(1)]="tumor"
# #   genomics$class.new[genomics$CNV.SNV_TP53 > 0]="tumor"
# #   table(genomics$class.new,useNA="ifany")
# #   # tumor    NA
# #   # 253     139
  
# #   #### . . . - Initializing a new object with the sample and class information - . . . ####
# #   summary.file.class=genomics[,c("sample","class.new")] ## reference file
  
# #   #### . . . - Defining tumor by counting the DRIVERs - . . . ####
# #   genomics_check=genomics[is.na(genomics$class.new),]
  
# #   for(i in 2:138){
# #     genomics_check[,i]=as.numeric(as.character(genomics_check[,i]))
# #     genomics_check[,i]=ifelse(is.na(genomics_check[,i]),0,genomics_check[,i])
# #   }
  
# #   #### . . . - Selecting the non-canonical drivers (without sampleID column) . . . - ####
# #   non.canonical.drivers=setdiff(colnames(genomics_check)[-c(1,139)],sig_feat)
  
# #   #### . . . - Creating all the possible variables to finish the classification . . . - ####
# #   genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-37]]) # exclude HY
# #   genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr21.gain")])# include HY and exclude all gains
# #   genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 37)]])# exclude HY and exclude all gains except for +18
# #   genomics_check$non_canonical_driver=rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains
# #   genomics_check$non_canonical_driver_no_gains=rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###
# #   genomics_check$only_gains=rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q
# #   genomics_check$igh=genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
# #     genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET
  
# #   #### . . . - Defining tumor: considering the presence of: if you have feature associated with SMM progression + hy or + igh translocation you are a tumor - . . . ####
# #   #dim(genomics_check) # 139
# #   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$igh > 0]="tumor" # IGH + at.least.1.driver
# #   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$hy > 0]="tumor"  # HY + at.least.1.driver 
# #   genomics_check$class.new[genomics_check$hy > 0 & genomics_check$igh > 0]="tumor"                      # HY + IGH
  
# #   #### . . . - Defining tumor: considering the presence of: if you have at least one MM prognostic features + any other feature you are tumor  - . . . ####
# #   # . - 1 is because the prognostic driver is counted twice (BAC: correct)
# #   genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$non_canonical_driver > 0]="tumor"
# #   #table(summary.file.class$class.new,useNA="ifany")
# #   summary.file.class$class.new[summary.file.class$sample %in% genomics_check$sample[genomics_check$class.new=="tumor"]]="tumor"
# #   #table(summary.file.class$class.new)
# #   # 297 tumor
  
# #   #### . . . - Defining non_tumor: considering the presence of: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
# #   # . - Define as non tumor all the samples with no drivers, but passed by filters
# #   genomics_check.2=genomics_check[is.na(genomics_check$class.new),]
# #   #dim(genomics_check.2)
# #   genomics_check.2$class.new[genomics_check.2$driver_count_all_gains == 0 &
# #                                genomics_check.2$driver_count_only_hy == 0 &
# #                                genomics_check.2$driver_count_only_no_hy == 0 &
# #                                genomics_check.2$non_canonical_driver == 0 &
# #                                genomics_check.2$non_canonical_driver_no_gains &
# #                                genomics_check.2$igh == 0]="non_tumor"
  
  
# #   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.2$class.new]="non_tumor"
# #   # table(summary.file.class$class.new,useNA="ifany")
# #   # tumor  <NA> 
# #   #   297    95
  
# #   #### . . . - Defining pts: check patients with more than one prognostic feature and define all the zeros pts as non tumor  - . . . ####
# #   genomics_check.3=genomics_check.2[is.na(genomics_check.2$class),]
# #   #dim(genomics_check.3) #95
# #   # genomics_check.3[genomics_check.3$driver_count_only_no_hy > 0,] ## all the samples with at least one genomic driver associated with progression are tumor
# #   genomics_check.3$non_canonical_driver_no_gains2=genomics_check.3$non_canonical_driver_no_gains - genomics_check.3$CNV.SNV_TGDS # we consider only one event in 13q
# #   #genomics_check.3[genomics_check.3$igh > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample # or >2
# #   #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sample #or >2
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
# #                                                                                                genomics_check.3$non_canonical_driver > 1))]="tumor"
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 1 |
# #                                                                                                 genomics_check.3$non_canonical_driver > 1))]="tumor"
# #   #genomics_check.3[genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0,]$sample
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$non_canonical_driver_no_gains2 > 1 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0)]="tumor"
  
# #   ##### presence of >2 non canonical is tumor - sse above
# #   #### presence of 2 + HY and IGH 
  
# #   #genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2),]$sample
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2))]="tumor" 
# #   #table(genomics_check.3$class.new,useNA="ifany")
# #   # tumor  <NA> 
# #   #   15    80
  
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains %in% c(0,1,2))]="non_tumor"  #c(0,1,2)
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(0,1,2))]="non_tumor" #c(0,1,2)
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0))]="non_tumor"
# #   genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0))]="non_tumor"
# #   #table(genomics_check.3$class.new,useNA="ifany")
# #   # non_tumor     tumor     
# #   #        80        15      
# #   #genomics_check.3[is.na(genomics_check.3$class.new),]
  
# #   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="tumor" ]]="tumor"
# #   summary.file.class$class.new[summary.file.class$sample %in% genomics_check.3$sample[genomics_check.3$class =="non_tumor" ]]="non_tumor"
  
# #   table(summary.file.class$class.new,useNA="ifany")
# #   # non_tumor     tumor 
# #   #        80       312 
  
# #   output$class=summary.file.class
# #   output$matrix=genomics[,-ncol(genomics)]
  
# #   return(output)
# # }


