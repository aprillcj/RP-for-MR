library("vroom")
library("data.table")
library("TwoSampleMR")
library("readxl")
library('gwasrapidd')
library("ggplot2")
library("cowplot")
library('dplyr')
library(biomaRt)
library(MungeSumstats)
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(gwasrapidd)
library(patchwork)
library(grid)
library(forestploter)
library(MRInstruments)

setwd("./Desktop/researchp/")
source("RP_func.R")


# pcos data:GCST007089.tsv DF
pcos<- read.csv("output.csv",header=F)
attributes(pcos)
pcos_outcome_dat<-format_data(pcos,type="outcome",
                    snp_col = "V2",
                    effect_allele_col = "V7",
                    other_allele_col = "V6",
                    eaf_col = "V8",
                    beta_col = "V9",
                    se_col = "V10",
                    pval_col = "V11",
                    chr_col = "V3",
                    phenotype = 'PCOS',
                    id = 'GCST007089'
)

pcos_exp_dat<-extract_instruments('finn-b-E4_POCS') # finn-b-E4-POCS
pcos_exp_dat2<-import_catalog("GCST007089.tsv") # GCST007089
pcos_exp_dat2$id.exposure='GCST007089'
pcos_exp_dat3<-extract_instruments('ebi-a-GCST90044902') # GCST90044902

# import the file from the excel file : ACCESSION ID 
sum <- read_xlsx("Summary_study.xlsx")
for(i in 1:nrow(sum)){
  loopid1<-sum$'GWAS_ID'[i]
  loopid1Nam<-sum$'Short_name'[i]
  loopid1Cat<-sum$'Catelogue'[i]
  print(paste("The loop for finn_4_PCOS:",loopid1,loopid1Nam,loopid1Cat))
  bi_ana_bat(loopid1,loopid1Nam,loopid1Cat)
  print(paste("The loop for bi_ana_bat2:",loopid1,loopid1Nam,loopid1Cat))
  bi_ana_bat2(loopid1,loopid1Nam,loopid1Cat)
  print(paste("The loop for bi_ana_bat3:",loopid1,loopid1Nam,loopid1Cat))
  bi_ana_bat3(loopid1,loopid1Nam,loopid1Cat)
}
bi_ana_bat3("ukb-b-17422","AAD2","Reproductive Trait")
bi_ana_bat3("ieu-b-4822","AAM","Reproductive Trait")
rm(exp_dat_AAM)
# PCOS: ebi-a-GCST90044902
# import the file from the ao file :ACCESSION ID 
ao <- available_outcomes() 
attributes(ao)
subao <- subset(ao,grepl("Barton AR",author,ignore.case = TRUE))
id_lsit <-subao$id
id_lsit<-id_lsit[1:3]
num=1
for (accessID in id_lsit) {
  print(paste0("Loop",num))
  if(exists(paste0("Res_tidy_",num))){
    print("ANALYSIS DO EXSITS")
  }
  else{
    print(accessID)
    print(paste("The loop for finn:",accessID,num))
    bi_ana_bat(accessID,num,"Barton AR-2021")
    print(paste("The loop for bi_ana_bat2:",accessID,num))
    bi_ana_bat2(accessID,num,"Barton AR-2021")
  }
  num <- num+1
}


########################################################
########################################################

# SUMRIZE THE DATA
all_objects <- ls()
objects <- grep("^Res_tidy", all_objects, value = TRUE) # for GWAS finn-b-E4_POCS
combined_res<-loop_append(objects) 



Method = list(
  "Inverse variance weighted" ="IVW",
  "MR Egger"="Egger",                                  
  "Weighted median"="WMedian",
  "Weighted mode"="WMode"
)

clarify<-function(x){
  if (grepl(': ', x)) { 
    gsub(".*: ", "", x)
  } else if (grepl("\\(.*\\)$", x)) {
    gsub("\\(.*", "", x)
  } else {
    x
  }
}
combined_res$item_method <- sapply(combined_res$exposure, function(x){
  if (grepl(': ', x)) { 
    gsub(".*: ", "", x)
  } else if (grepl("\\(.*\\)$", x)) {
    gsub("\\(.*", "", x)
  } else {
    x
  }
})


library("writexl")
write_xlsx(combined_res,"op/result3.xlsx")
write_xlsx(combined_test,"op/result_test0.xlsx")

# PCOS AS EXPOSURE
# combined_res$item_method <- paste(combined_res$exposure,"\n",Method[combined_res$method])

#finn-b-E4_POCS
combined_res<-subset(combined_res,exposure!='Waist-hip ratio')

subset_BMI = combined_res[combined_res$id.exposure %in% c('ieu-a-974', 'ieu-a-835', 'ieu-a-785'), ]
subset_not_BMI = combined_res[!combined_res$id.exposure %in% c('ieu-a-974', 'ieu-a-835', 'ieu-a-785'),  ]

subset_BMI$item_method <- ifelse(subset_BMI$id.exposure == 'ieu-a-974', paste(subset_BMI$item_method, "total"), subset_BMI$item_method)
subset_BMI$item_method <- ifelse(subset_BMI$id.exposure == 'ieu-a-835', paste(subset_BMI$item_method, "female"), subset_BMI$item_method)
subset_BMI$item_method <- ifelse(subset_BMI$id.exposure == 'ieu-a-785', paste(subset_BMI$item_method, "male"), subset_BMI$item_method)

combined_IVW<-subset(subset_not_BMI,method=="Inverse variance weighted"&id.outcome=="finn-b-E4_POCS" )
Visualizaton(combined_IVW,"PCOS:finn-b-E4_POCS AS OUTCOME")

sss_BMI<-subset(subset_BMI,method=="Inverse variance weighted"&id.outcome=="finn-b-E4_POCS" )
visual_simple(sss_BMI,"BMI")

# ebi-a-GCST90044902

combined_IVW2<-subset(subset_not_BMI,method=="Inverse variance weighted"&id.outcome=="ebi-a-GCST90044902" )
Visualizaton(combined_IVW2,"PCOS:ebi-a-GCST90044902 AS OUTCOME")
t_combined_IVW2<-subset(subset_not_BMI,id.outcome=="ebi-a-GCST90044902" )

Pos_result$comment <- with(Pos_result, paste0(item_method, 
                                             ' (', 
                                             round(or, 2), 
                                             '[', 
                                             round(or_lci95, 2), 
                                             "-", 
                                             round(or_uci95, 2), 
                                             ']; P=', 
                                             format.pval(pval, digits = 3),
                          
                                             ')'))
Pos_result$comment
notsig<-subset(combined_IVW2,pval>0.05)
notsig$comment <- with(notsig, paste0(item_method, 
                                      ' (', 
                                      'P=', 
                                      format.pval(pval, digits = 3),
                                      ')'))
sorted_nonsig_df <- notsig[order(notsig$pval), ]
sorted_nonsig_df$comment


sss_BMI2<-subset(subset_BMI,method=="Inverse variance weighted"&id.outcome=="ebi-a-GCST90044902" )
visual_simple(sss_BMI2,"BMI2")



# PCOS as the exposure
object2<- grep("^Bi_Res_tidy",  ls(), value = TRUE) 
combined_bi_res<-loop_append(object2)


combined_bi_res$item_method <- sapply(combined_bi_res$outcome, function(x){
  if (grepl(': ', x)) { 
    gsub(".*: ", "", x)
  } else if (grepl("\\(.*\\)$", x)) {
    gsub("\\(.*", "", x)
  } else {
    x
  }
})

subset_BMI2 = combined_bi_res[combined_bi_res$id.outcome %in% c('ieu-a-974', 'ieu-a-835', 'ieu-a-785'), ]
subset_not_BMI2 = combined_bi_res[!combined_bi_res$id.outcome %in% c('ieu-a-974', 'ieu-a-835', 'ieu-a-785'),  ]

combined_bi_IVW<-subset(subset_not_BMI2,method=="Inverse variance weighted"&id.exposure=="ebi-a-GCST90044902")
Visualizaton(combined_bi_IVW,"PCOS:ebi-a-GCST90044902 AS EXPOSURE")

Visualizaton(combined_IVW2,"PCOS AS EXPOSURE")

library(MRMix)

Pos_result <- subset(combined_IVW2,pval <= 0.05)
Pos_result$name<- sapply(Pos_result$resource,function(x){ gsub(".*_", "", x)})
Pos_result$type <- ifelse(Pos_result$item_method == "Age at menopause", "binary", "continuous")

for(i in 1:nrow(Pos_result)){
  id<-Pos_result$'id.exposure'[i]
  idname<-Pos_result$'name'[i]
  idtype<-Pos_result$'type'[i]
  print(paste(id,idname))
  newanalysis(id,idname,idtype)
}

objects <- grep("^Check", ls(), value = TRUE) 
combined_res<-loop_append(objects) 



# combined senetivity check
objects <- grep("^Out_sen3", ls(), value = TRUE) 
combined_sen<-loop_append(objects) 

gsub(".*||", "", x)

combined_sen$item_method <- sapply(combined_sen$exposure, function(x){
  y<-gsub("\\|\\| .*", "", x)
})
combined_sen$id.exposure <- sapply(combined_sen$exposure, function(x){
  gsub(".*id\\:", "", x)
})

op_sen<-combined_sen[Pos_result$id.exposure]
op_sen <- combined_sen[combined_sen$id.exposure %in% Pos_result$id.exposure, ]
write_xlsx(op_sen,"op/sen3.xlsx")

