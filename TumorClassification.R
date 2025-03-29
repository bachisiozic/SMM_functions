



data=read.delim("~/Desktop/BACHI/Project1/00_FINALdata/2025_FINAL_03172025.txt",stringsAsFactors = F)
data[data$sampleID == "52915_DNA_T",]

sig_feat=c("CNV_chr3.gain","CNV_chr5.gain","CNV_chr7.gain","CNV_chr9.gain","CNV_chr11.gain","CNV_chr15.gain","CNV_chr19.gain","CNV_chr18.gain",
           "CNV_Gain_Amp1q",
           "CNV_chr8p.loss",
           "CNV_Del_10q24.32","CNV_Del_10q26.3","CNV_Del_12q24.31","CNV_Del_14q24.3","CNV_Del_2q31.1","CNV_Del_2q37.3","CNV_Del_6q26","CNV_Del_8q24.21",
           "CNV.SNV_ARID2","CNV.SNV_BTG1","CNV.SNV_DNMT3A","CNV.SNV_TENT5C","CNV.SNV_NCOR1","CNV.SNV_NF1","CNV.SNV_POT1","CNV.SNV_PRDM1","CNV.SNV_TET2",
           "CNV.SNV_TP53","CNV.SNV_ZNF292",
           "SNV_FGFR3","SNV_NRAS","SNV_BHLHE41","SNV_HIST1H1D",
           "APOBEC",
           "t_MMSET",
           "CNV.Sig",
           "hy")

dim(data)
#392
str(data)

copynumbers=colnames(data)[grep("CNV_",colnames(data))]
tumorsuppressorgenes=colnames(data)[grep("CNV.SNV_",colnames(data))]
oncogenesothers=c(colnames(data)[grep("^SNV_|^SNV.",colnames(data))])
genomics=data[,c(colnames(data)[grep("sample|sampleID|Sample|SampleID|SAMPLE|SAMPLE.ID",colnames(data))],
                 copynumbers,tumorsuppressorgenes,oncogenesothers,
                 "CNV.Sig","APOBEC","t_CCND1","t_MMSET","t_MAF","CCND3","CCND2","MAFA","MYC","hy")]



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

genomics[genomics$sampleID == "Cap1443-4-ID09",]
# FINAL$class[FINAL$sampleID == "NIH032A"]<- "tumor" ### qc MOSTRA COMPLEX on 16p and 16q AND 11;14
# FINAL$class[FINAL$sampleID == "78208_DNA_T"]<- "tumor"
# FINAL[FINAL$sampleID =="NIH036A",] # check APOBEC right convertion
# FINAL[FINAL$sampleID =="BM046361",]
# FINAL$class[FINAL$sampleID == "Cap1443-4-ID09"]<-"tumor"

######################################
##
#### . . . - DEFINING TUMOR BY COUNTING THE DRIVERs
##
######################################

summary.file.class=genomics[,c("sampleID","class.new")] ## reference file

genomics_check=genomics[is.na(genomics$class.new),]

for(i in 2:138){
  genomics_check[,i]=as.numeric(as.character(genomics_check[,i]))
  genomics_check[,i]=ifelse(is.na(genomics_check[,i]),0,genomics_check[,i])
}


# genomics_check[,2:137]=apply(genomics_check[,2:137], 2, function(x){as.numeric(as.character(x))})
# genomics_check<- genomics_check[!is.na(genomics_check$sampleID),]
# genomics_check$APOBEC[is.na(genomics_check$APOBEC)]<-0
# genomics_check$CCND2[is.na(genomics_check$CCND2)]<-0
# genomics_check$CCND3[is.na(genomics_check$CCND3)]<-0
# genomics_check$MAFA[is.na(genomics_check$MAFA)]<-0
# genomics_check$CCND2[is.na(genomics_check$t_CCND1)]<-0
# genomics_check$CCND3[is.na(genomics_check$t_MMSET)]<-0
# genomics_check$MAFA[is.na(genomics_check$t_MMSET)]<-0
# sig_feat2<- sig_feat[-which(sig_feat =="Simple")] ## remove Simple


#### . . . - Selecting the non-canonical drivers (without sampleID column)
non.canonical.drivers=setdiff(colnames(genomics_check)[-c(1,139)],sig_feat)


genomics_check$driver_count_all_gains=rowSums(genomics_check[,sig_feat[-37]]) # exclude HY
genomics_check$driver_count_only_hy=rowSums(genomics_check[,c(sig_feat[c(1:7)],"CNV_chr21.gain")])# include HY and exclude all gains
genomics_check$driver_count_only_no_hy=rowSums(genomics_check[,sig_feat[-c(1:7, 37)]])# exclude HY and exclude all gains except for +18

genomics_check$non_canonical_driver<- rowSums(genomics_check[,c(non.canonical.drivers[-c(1,95:100)],"CNV_chr18.gain")]) # chr 18 included as well as all non HY gains

genomics_check$non_canonical_driver_no_gains<- rowSums(genomics_check[,c(non.canonical.drivers[c(8:94)])]) ###chr##18##included###

genomics_check$only_gains<- rowSums(genomics_check[,c(2:10, 12:17)]) # all possible gains except 1q

genomics_check$igh<- genomics_check$MAFA + genomics_check$CCND2 + genomics_check$CCND3 + 
  genomics_check$t_CCND1 + genomics_check$t_MAF + genomics_check$t_MMSET



########################################################################################################################
#
# if you have feature associated with SMM progression + hy or + igh translocation you are a tumor
###
###########################################################################################################################
table(genomics_check$class.new)
dim(genomics_check) # 139

genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$igh > 0]="tumor" # IGH + at.least.1.driver
genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$hy > 0]="tumor"  # HY + at.least.1.driver 
genomics_check$class.new[genomics_check$hy > 0 & genomics_check$igh > 0]="tumor"                      # HY + IGH

# # colnames(FINAL)[which(as.numeric(FINAL[FINAL$sampleID =="SL321965",10:ncol(FINAL)])>0)+9]
# 
# as.character(na.omit(FINAL_check$sampleID[FINAL_check$disease_stage=="MGUS" & FINAL_check$class=="tumor"]))
# # "SL235706"
# 
# table(FINAL_check$disease_stage, FINAL_check$class)

############################################################################################################
###
## 2. if you have at least one MM prognostic features + any other feature you are tumor 
# 1 is because the prognostic driver is counted twice (BAC: correct)
###
##################################################################################################################
genomics_check$class.new[genomics_check$driver_count_only_no_hy > 0 & genomics_check$non_canonical_driver > 0]="tumor"
# FINAL_check[FINAL_check$sampleID == "52915_DNA_T",]

table(summary.file.class$class.new,useNA="ifany")
summary.file.class$class.new[summary.file.class$sampleID %in% genomics_check$sampleID[genomics_check$class.new=="tumor"]]="tumor"
table(summary.file.class$class.new,useNA="ifany")


############################################################################################################
###
### 3. check patients with more than one prognostic feature and define all the zeros pts as non tumor
###
##################################################################################################################
genomics_check.2=genomics_check[is.na(genomics_check$class.new),]
dim(genomics_check.2)
# genomics_check.2[genomics_check.2$driver_count_only_no_hy > 0 ,c(139:ncol(genomics_check.2))]
# # "SL235769" - this is a non tumor
# # "SL332283" -  this is a low purity - 5% VAF likely high purity - this sample has been removed above. 
# FINAL_check2$class[FINAL_check2$sampleID == "SL235769"]<-"non_tumor"
# # FINAL_check2<- FINAL_check2[(FINAL_check2$sampleID) != "SL332283",] # remove bad sample - it is remove above as well
# 
# summary_file_class$class[summary_file_class$sampleID == "SL235769"]<-"non_tumor"


## define as non tumor all the samples with no drivers, but passed by bachi and me
genomics_check.2$class.new[genomics_check.2$driver_count_all_gains == 0 &
                             genomics_check.2$driver_count_only_hy == 0 &
                             genomics_check.2$driver_count_only_no_hy == 0 &
                             genomics_check.2$non_canonical_driver == 0 &
                             genomics_check.2$non_canonical_driver_no_gains &
                             genomics_check.2$igh == 0]="non_tumor"


summary.file.class$class.new[summary.file.class$sampleID %in% genomics_check.2$class.new]="non_tumor"
table(summary.file.class$class.new,useNA="ifany")


############################################################################################################
###
### 4. check patients with more than one prognostic feature and define all the zeros pts as non tumor
###
##################################################################################################################
genomics_check.3=genomics_check.2[is.na(genomics_check.2$class),]
dim(genomics_check.3) #95

genomics_check.3[genomics_check.3$driver_count_only_no_hy > 0,] ## all the samples with at least one genomic driver associated with progression are tumor at this point

genomics_check.3$non_canonical_driver_no_gains2=genomics_check.3$non_canonical_driver_no_gains - genomics_check.3$CNV.SNV_TGDS # we consider only one event in 13q

genomics_check.3[genomics_check.3$igh > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sampleID # or >2
#### BACHISIO
# "BM041752"
# "Cap1443-4-ID09"
# "NIH010A"
# "SL217777"  

genomics_check.3[genomics_check.3$sampleID == "SL217777",] # tumor
genomics_check.3[genomics_check.3$sampleID == "BM041752",] # tumor
genomics_check.3[genomics_check.3$sampleID == "Cap1443-4-ID09",] # tumor
genomics_check.3[genomics_check.3$sampleID == "NIH010A",] # tumor

#### FRANCESCO
# "BM041752" -->"CNV_Del_10p15.3" "CNV_Del_7p22.2"  "CNV.SNV_FUBP1"   "t_CCND1"         "Simple"  --> ???? tumoral 
# "SL217777" --> "CNV.SNV_RB1"  "CNV.SNV_TGDS" "SNV_CCND1"    "t_CCND1"      "Simple"  --> TUMOR!

genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 > 1,]$sampleID #or >2
#### BACHISIO
# "228503-WX01"
# "Cap1643-4-ID24"
# "SL240420"
# "SL241594"
# "SRR18145704"
# "SRR18145707"
#### FRANCESCO
# "228503-WX01" 
# "SL240420"    
# "SRR18145707"

genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 2 |
                                                                                             genomics_check.3$non_canonical_driver > 2))]="tumor"
genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 & (genomics_check.3$non_canonical_driver_no_gains2 > 2 |
                                                                                              genomics_check.3$non_canonical_driver > 2))]="tumor"

genomics_check.3[genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0,]$sampleID
#### BACHISIO
# "228509-WX01"
# "PD57074a"
# "SRR18145694"
# "SRR18145696"
#### FRANCESCO
# "228509-WX01"
# "228519-WX01"
# "PD57074a" - 20q del and 13q and cosmic genes

genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$non_canonical_driver_no_gains2 > 2 & genomics_check.3$hy == 0 & genomics_check.3$igh == 0)]="tumor"

##### presence of >2 non canonical is tumor - sse above

#### presence of 2 + HY and IGH 

genomics_check.3[genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2),]$sampleID
#### BACHISIO 
# "Cap1643-4-ID24" 
# "SL241594"
# "SRR18145704"  

genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2))]="tumor" 

genomics_check.3[genomics_check.3$igh > 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(2),]$sampleID
# "NIH010A" 
# "SL217777"

table(genomics_check.3$class.new,useNA="ifany")
# tumor  <NA> 
#   18    77

genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains %in% c(0,1,2))]="non_tumor"  #c(0,1,2)
genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy == 0 & genomics_check.3$igh == 0 & genomics_check.3$non_canonical_driver_no_gains2 %in% c(0,1,2))]="non_tumor" #c(0,1,2)
genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$hy > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
genomics_check.3$class.new[is.na(genomics_check.3$class.new) & (genomics_check.3$igh > 0 &  genomics_check.3$non_canonical_driver_no_gains2 %in% c(1,0,2))]="non_tumor"
table(genomics_check.3$class.new,useNA="ifany")
# non_tumor     tumor     
#        77        18      

genomics_check.3[is.na(genomics_check.3$class.new),]

# genomics_check.3$sampleID[is.na(genomics_check.3$class)]
# # character(0)
# 
# table(genomics_check.3$class.new)## none of these patient progress
# 
# ## nonly two pts imwg 2
# table(FINAL_check3$imwg)
# FINAL_check3[FINAL_check3$imwg %in% c("3", "2"),]
# 
# # Cap1453-4-ID19 - IMWG2 with 2670 days no progression
# 
# # 226759-AU01 low purity --> to be removed
# #

summary.file.class$class.new[summary.file.class$sampleID %in% genomics_check.3$sampleID[genomics_check.3$class =="tumor" ]]="tumor"
summary.file.class$class.new[summary.file.class$sampleID %in% genomics_check.3$sampleID[genomics_check.3$class =="non_tumor" ]]="non_tumor"
# summary_file_class[is.na(summary_file_class$class)]
# data frame with 0 columns and 395 rows
summary.file.class[summary.file.class$sampleID == "Cap1453-4-ID19",]
summary.file.class[summary.file.class$sampleID == "NIH010A",]

table(summary.file.class$class.new,useNA="ifany")
# non_tumor     tumor 
#        77       315 
#summary.file.class$class[summary.file.class$sampleID == "NIH010A"]="non_tumor"


# table(summary_file_class$disease_stage, summary_file_class$class)
# #          non_tumor tumor
# # MGUS        46    43
# # SMM         29   274
# # 
write.table(summary.file.class,
            "~/Desktop/BACHI/Project1/GitHub/SMM/2025_tumor_vs_nontumor.txt",
            col.names = T,
            row.names = F,
            quote = F,
            sep="\t")



########### 
colnames(summary.file.class)[2]="class.new.BAC"

data.new=read.delim("~/Desktop/BACHI/Project1/00_FINALdata/2025_tumor_vs_nontumor_03172025.txt",stringsAsFactors = F)
dim(data.new) # 392
dim(summary.file.class) # 392

def=merge(data.new,summary.file.class,by="sampleID")
dim(def)

table(def$class,useNA="ifany")
table(def$class.new.BAC,useNA="ifany")
table(def$class,def$class.new.BAC,useNA="ifany")

def$FILTER=ifelse(def$class == def$class.new.BAC,"MATCH","NO")
table(def$FILTER)
# MATCH    NO 
#   384     5 

def[def$FILTER=="NO",]
#           sampleID disease_stage class class.new.BAC FILTER
# 115 Cap1453-4-ID19           SMM tumor     non_tumor     NO # maybe to remove it (From the QC it seems good)
# 128 Cap2263-4-ID09          MGUS tumor     non_tumor     NO
# 129 Cap2264-4-ID10           SMM tumor     non_tumor     NO
# 329       SL321965          MGUS tumor     non_tumor     NO
# 343       SL380854          MGUS tumor     non_tumor     NO

genomics_check.3[genomics_check.3$sampleID == "Cap1453-4-ID19",]
genomics_check.3[genomics_check.3$sampleID == "Cap2263-4-ID09",]
genomics_check.3[genomics_check.3$sampleID == "Cap2264-4-ID10",]
genomics_check.3[genomics_check.3$sampleID == "SL321965",]
genomics_check.3[genomics_check.3$sampleID == "SL380854",]

genomics_check.3[genomics_check.3$sampleID %in% def[def$FILTER=="NO","sampleID"],c("sampleID","hy","class.new")]
fra2[fra2$sampleID %in% def[def$FILTER=="NO","sampleID"],c("sampleID","hy","pfs_code","pfs_time")]


table(genomics_check.3$hy,genomics_check.3$class.new)





