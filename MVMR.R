


# telestero
id_exposure <- c("ebi-a-GCST90012114","ebi-a-GCST90012104")
id_outcome <- "finn-b-E4_POCS" 
exposure_dat<- mv_extract_exposures(id_exposure) 
oucome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
mvdat <- mv_harmonise_data(exposure_dat, oucome_dat)  
mvres <- mv_multiple(mvdat) 


id_exposure <- c("ieu-b-4824","finn-b-N14_MESNRUIRREG")
exposure_dat2<- mv_extract_exposures(id_exposure)
oucome_dat2 <- extract_outcome_data(exposure_dat2$SNP, id_outcome)
mvdat_2 <- mv_harmonise_data(exposure_dat2, oucome_dat2)  
mvres2 <- mv_multiple(mvdat_2) 

#Glycemic
id_exposure <- c("ebi-a-GCST90002238","ebi-a-GCST90002232")
id_exposure <- c("ebi-a-GCST90014006","ieu-b-117")
exposure_dat3<- mv_extract_exposures(id_exposure) 
oucome_dat3 <- extract_outcome_data(exposure_dat3$SNP, id_outcome)
mvdat3 <- mv_harmonise_data(exposure_dat3, oucome_dat3)  
mvres3 <- mv_multiple(mvdat3) 
res_Hba1c$dataset2
res_HOMAB$dataset2

id_exposure <- c("ebi-a-GCST90013974","ieu-b-110")
exposure_dat4<- mv_extract_exposures(id_exposure)
oucome_dat4 <- extract_outcome_data(exposure_dat4$SNP, id_outcome)
mvdat_4 <- mv_harmonise_data(exposure_dat4, oucome_dat4)  
mvres4 <- mv_multiple(mvdat_4) 
res_LDL$dataset1


# IGF-1 and BMI
id_exposure <- c("ebi-a-GCST90025989","ebi-a-GCST90025994")
exposure_dat5<- mv_extract_exposures(id_exposure)
oucome_dat5 <- extract_outcome_data(exposure_dat5$SNP, id_outcome)
mvdat_5 <- mv_harmonise_data(exposure_dat5, oucome_dat5)  
mvres5 <- mv_multiple(mvdat_5) 
?mv_multiple
aobmi<-subset(ao,grepl("body mass index",trait,ignore.case = T))

# Total bilirubin levels and BMI 
id_exposure<-c("ebi-a-GCST90014012","ieu-b-4815")
exposure_dat6<- mv_extract_exposures(id_exposure)
oucome_dat6 <- extract_outcome_data(exposure_dat6$SNP, id_outcome)
mvdat_6 <- mv_harmonise_data(exposure_dat6, oucome_dat6)  
mvres6 <- mv_multiple(mvdat_6) 
mvres6 <- mv_ivw(mvdat_6) 

# Serum albumin levels and BMI 
id_exposure <- c("ebi-a-GCST90025992","ebi-a-GCST90025973","ieu-b-4815")
exposure_dat7<- mv_extract_exposures(id_exposure)
oucome_dat7 <- extract_outcome_data(exposure_dat7$SNP, id_outcome)
mvdat_7 <- mv_harmonise_data(exposure_dat7, oucome_dat7)  
mvres7 <- mv_multiple(mvdat_7) 
mvres7<-mv_ivw(mvres7)

#Bilirubin levels , Direct bilirubin levels , Body mass index 
id_exposure <- c("ebi-a-GCST90025983","ebi-a-GCST90025973","ebi-a-GCST90025994")
exposure_dat8<- mv_extract_exposures(id_exposure)
oucome_dat8 <- extract_outcome_data(exposure_dat8$SNP, id_outcome)
mvdat_8 <- mv_harmonise_data(exposure_dat8, oucome_dat8)  
mvres8 <- mv_multiple(mvdat_8) 
mvres8<-mv_ivw(mvdat_8)
#
#Bilirubin levels , Direct bilirubin levels 
id_exposure <- c("ebi-a-GCST90025983","ebi-a-GCST90025973")
exposure_dat9<- mv_extract_exposures(id_exposure)
oucome_dat9 <- extract_outcome_data(exposure_dat9$SNP, id_outcome)
mvdat_9 <- mv_harmonise_data(exposure_dat9, oucome_dat9)  
mvres9 <- mv_multiple(mvdat_9) 
mvres9<-mv_ivw(mvdat_9)
mvres9$plots[1]
mvres9$plots[2]

id_exposure <- c("ieu-a-835","ieu-a-785")
exposure_dat9<- mv_extract_exposures(id_exposure)
oucome_dat9 <- extract_outcome_data(exposure_dat9$SNP, id_outcome)
mvdat_9 <- mv_harmonise_data(exposure_dat9, oucome_dat9)  
mvres9 <- mv_multiple(mvdat_9) 
mvres9<-mv_ivw(mvdat_9)


id_exposure <- c("ieu-b-4864","prot-a-529","ebi-a-GCST004941","ebi-a-GCST005179") 
id_exposure <- c("ieu-b-4864","prot-a-529","ebi-a-GCST004941") 

exposure_dat <- mv_extract_exposures(id_exposure) 
id_exposure <- c("ieu-b-4816","ieu-b-102") # depression & bmi
exposure_dat <- mv_extract_exposures(id_exposure)  
print(dim(exposure_dat))
# AMH_exp_dat <-obtain_gwas(pheno_key="anti-mullerian")
# col_list<- colnames(exposure_dat)
# AMH_exp_dat<-AMH_exp_dat[,col_list]
# AMH_exp_dat$exposure<-paste0(AMH_exp_dat$exposure,'|| id:',AMH_exp_dat$id.exposure)
# keepsnps <- exposure_dat$SNP
# AMH_exp_dat2 <- subset(AMH_exp_dat, SNP %in% keepsnps)
# exposure_dat1<-rbind(exposure_dat,AMH_exp_dat2)

col_list<- colnames(oucome_dat)
colnames(pcos_o_dat)
pcos_o_dat<-pcos_o_dat[,col_list]

test <- extract_instruments("finn-b-E4_POCS") 
dim(test)
oucome_dat <- extract_outcome_data(exposure_dat$SNP, 'finn-b-E4_POCS')

mvdat <- mv_harmonise_data(exposure_dat, oucome_dat) 

dim(exposure_dat)
dim(oucome_dat)
attributes(mvdat)
# debug of the mv_multiple
beta.outcome <- mvdat$outcome_beta
beta.exposure <- mvdat$exposure_beta
w <- 1/mvdat$outcome_se^2
mod <- summary(stats::lm(beta.outcome ~ 0 + 
                           beta.exposure, weights = w))
mod$coef


res1 <- mv_multiple(mvdat)
res <- mv_multiple(mvdat,instrument_specific=TRUE)
OR <- generate_odds_ratios(res1) 

mvdat1 <- mv_harmonise_data(exposure_dat, oucome_dat)
attributes(mvdat1)
beta.outcome <- mvdat1$outcome_beta
beta.exposure <- mvdat1$exposure_beta
w <- 1/mvdat1$outcome_se^2
mod <- summary(stats::lm(beta.outcome ~ 0 + 
                           beta.exposure, weights = w))
mod$coef
res1 <- mv_multiple(mvdat1)

mvres$result
res_T_total$dataset1


res_BMI$dataset2


id_exposure <- c("ebi-a-GCST90012114","ebi-a-GCST90012104")
id_exposure <- c("ieu-b-109","ebi-a-GCST90013974")
id_outcome <-"ebi-a-GCST90044902"
exposure_dat1<- mv_extract_exposures(id_exposure)
oucome_dat1 <- extract_outcome_data(exposure_dat1$SNP, id_outcome)
mvdat <- mv_harmonise_data(exposure_dat1, oucome_dat1)  
mvres9 <- mv_multiple(mvdat) 
mvres9<-mv_ivw(mvdat)

ieu-b-109 HDL
 
chd_out_dat <- extract_outcome_data(
    snps = exp_dat_HDL$SNP,
    outcomes = "ebi-a-GCST90013974")

dat <- harmonise_data(
    exposure_dat = exp_dat_HDL, 
    outcome_dat = chd_out_dat
  )
test<-mr(dat)

chd_out_dat <- extract_outcome_data(
  snps = exp_dat_BMI$SNP,
  outcomes = "ieu-b-109")

dat <- harmonise_data(
  exposure_dat = exp_dat_BMI, 
  outcome_dat = chd_out_dat
)
test<-mr(dat)
mvres1
comeid<-"ebi-a-GCST90044902"
idlist2 <- c("ebi-a-GCST90002238","ebi-a-GCST90013974")
mv2<-multiple_ana(idlist2,comeid)

ukb-b-17422
comeid<-"ebi-a-GCST90044902"
idlist2 <- c("ebi-a-GCST90018959","ebi-a-GCST90013974")
mv4<-multiple_ana(idlist2,comeid)



