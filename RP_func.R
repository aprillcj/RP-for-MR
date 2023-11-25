
#ANALYSIS FOR BIDIRECTION
bi_ana_bat<-function(id1,id1Nam,id1Cat){
  # id1<-"ebi-a-GCST90013974"
  # id1Nam<-38
  # id1Cat<-"Physical Trait"
  ###########################################PCOS AS OUTCOME#########################################
  # exact the exposure from the IEU DATABASE
  print("START ANALYSIS PCOS AS OUTCOME")
  name <- paste0('exp_dat_',id1Nam)
  assign(eval(name),extract_instruments(eval(id1)),envir = .GlobalEnv)
  
  if(!exists(eval(name))){
    print("ACHIEVE ASSOCIATION FAILED")
    return("NO ASSOCIATION SNPs")
  }else{
    print("expDS")
    expDS <- get(eval(name))
  }
  
  print('chd_our_dat')

  chd_out_dat <- tryCatch(
    extract_outcome_data(
      snps = expDS$SNP,
      outcomes = "finn-b-E4_POCS"
    ),
    error=function(e){
      # print("chd_out_dat ERROR")
      stop("chd_out_dat ERROR")
    }
  )
  
  if (is.null(chd_out_dat)) {
    cat("Skipping analysis due to missing chd_out_dat.\n")
    return(NULL)
  }
  
  
  print("harmonise")
  # harmonise data
  dat <- tryCatch(
    harmonise_data(
      exposure_dat = expDS, 
      outcome_dat = chd_out_dat
    )
    ,
    error=function(e){
      cat("Error in harmonise_data. Skipping to the next iteration.\n")
      return(NULL)
      
    }
  )
  if (is.null(dat)) {
    cat("Skipping analysis due to missing data.\n")
    return(NULL)
  }
  
  
  #analysis 
  res<-mr(dat)
  print('res finish')
  
  if(nrow(res)==0){
    print("RES 0 ROW")
    return("RES 0 ROW")
  }
  
  # Tidy up results and save
  res_tidy<-res %>%
    split_outcome() %>% 
    split_exposure() %>%
    generate_odds_ratios()
  print("analysis finish")
  # return(res_tidy)
  res_tidy$'factor'<-id1Cat
  assign(paste0("Res_tidy_",id1Nam), res_tidy, envir = .GlobalEnv)
  
  if(nrow(dat)<=3){
    return("dat rows <= 3 ")
  }
  
  print("Heterogeneity and Pleiotropy")
  mr_heterogeneity(dat,method_list=c("mr_ivw_mre","mr_egger_regression"))
  out_sensitivity <- full_join(mr_pleiotropy_test(dat), 
                               mr_heterogeneity(dat)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list <- list(dataset1 = res_tidy, dataset2 = out_sensitivity)
  assign(paste0("Out_sen_",id1Nam), out_sensitivity, envir = .GlobalEnv)
  
  print("Sensitivity analysis")
  # sensitivity analysis
  sen <- mr_leaveoneout(dat)
  
  
  ###########################################PCOS AS EXPOSURE#########################################
  print("START ANALYSIS PCOS AS EXPOSURE")
  
  chd_out_dat2 <- tryCatch(
    extract_outcome_data(
      snps = pcos_exp_dat$SNP,
      outcomes = id1
    ),
    error=function(e){
      print("chd_out_dat ERROR")
      # stop("chd_out_dat ERROR")
    }
  )
  
  if (is.null(chd_out_dat2)) {
    cat("Skipping analysis due to missing chd_out_dat.\n")
    return(NULL)
  }
  
  dat2 <- harmonise_data(
    exposure_dat = pcos_exp_dat,
    outcome_dat = chd_out_dat2
  )
  if (is.null(dat2)) {
    cat("Skipping analysis due to dat2.\n")
    return(NULL)
  }
  
  res2<-mr(dat2)
  res_tidy2<-res2 %>%
    split_outcome() %>%
    split_exposure() %>%
    generate_odds_ratios()
  res_tidy2$'factor'<-id1Cat
  assign(paste0("Bi_Res_tidy_",id1Nam), res_tidy2, envir = .GlobalEnv)
  
  
  
  if(nrow(dat2)<=3){
    return("dat2 rows <= 3 ")
  }
  out_sensitivity2 <- full_join(mr_pleiotropy_test(dat2),
                                mr_heterogeneity(dat2)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list2 <- list(dataset1 = res_tidy2, dataset2 = out_sensitivity2)
  assign(paste0("Bi_Out_sen_",id1Nam), out_sensitivity2, envir = .GlobalEnv)
  
  print('FINISH BIDIRECTION ANALYSIS')
  
}


# Visualizaton 
# Import the forest plot
Visualizaton<-function(Obj,name){
  pal<-c(unname(yarrr::piratepal("pony")))
  Obj$pcate <- ifelse(Obj$pval < 0.05, "<0.05",
                      ifelse(Obj$pval <0.1, "0.05-0.1", ">=0.1"))
  min<-min(Obj$or_lci95)-0.5
  max<-max(Obj$or_uci95)+0.5
  Obj <- Obj %>%
    arrange(factor)

  p1<-ggplot(Obj, aes(y=item_method,x=or),size=4,title="xx",label="te") +
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95,col=factor),height=.3) +
    geom_point(aes(size=se,shape= pcate,col=factor))+
    scale_color_manual(values=pal)+
    # xlim(min,max)+
    geom_vline(xintercept=1, linetype='longdash',color="blue") +
    theme_minimal_grid(9) +
    facet_grid(factor~.,scale="free",switch = "y",space="free") +
    geom_text(aes(x = max+1,label = sprintf("%0.2f", round(or, digits = 2)))) +
    geom_text(aes(x = max+2,label = sprintf("%0.2f", round(b, digits = 2)))) + 
    geom_text(aes(x = max+3,label = sprintf("%1.f", round(nsnp, digits = 2)))) + 
    # scale_y_discrete(limits = rev(levels(Obj$item_method))) +
    scale_x_continuous(
      breaks =  (max+1):(max+3) , 
      labels = c('OR', 'Beta','nSNP'), 
      position = 'top' ) +
    theme(strip.text = element_text(size=4,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom" ,
          # legend.text = element_text(size=5),
          plot.title = element_text(hjust = 0.5)
          
    )+
    guides(
      col = guide_legend(title = "Factor",ncol=1),
      size = guide_legend(title = "SE Size",ncol=1),
      shape = guide_legend(title = "Pvalue",ncol=1)
    )+
    labs(title = "",x="",y="",subtitle = "")
    
  
  
  p2<-ggplot(Obj, aes(y = item_method)) +
    geom_text(aes(x = 0, label = sprintf("%0.3f", round(or, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = sprintf("%0.3f", round(b, digits = 2))), size = 4) + 
    geom_text(aes(x = 2, label = sprintf("%1.f", round(nsnp, 1))), size = 4) + 
    # scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_y_discrete(limits = rev(levels(Obj$item_method))) +
    scale_x_continuous(
      breaks = 0:2, labels = c('OR', 'Beta','nSNP'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))

  p <- p1 + p2 + plot_layout(widths = c(0.78, 0.2))
  
  total_height <- nrow(Obj) * 0.3
  ggsave(paste0("op/", name,"_forest.pdf"),
         plot = p1,
         width = 8,
         height = total_height, 
         units = "in",
         dpi = 600,
         limitsize = FALSE)
  ggsave(paste0("op/",name, "_forest.jpg"),
         plot = p1,
         width = 8,
         height = total_height,
         units = "in",
         dpi = 600,
         limitsize = FALSE)
  
  print('The picture has been saved')
  
}



bi_ana_bat2<-function(id1,id1Nam,id1Cat){
  # id1<-"ebi-a-GCST90013974"
  # id1Nam<-38
  # id1Cat<-"Physical Trait"
  ###########################################PCOS AS OUTCOME#########################################
  # exact the exposure from the IEU DATABASE
  print("START ANALYSIS PCOS AS OUTCOME")

  name <- paste0('exp_dat_',id1Nam)
  if(exists(name)){
    print(paste(name,"has been exists"))
  }else{
    assign(eval(name),extract_instruments(eval(id1)),envir = .GlobalEnv)
  }
  
  if(!exists(eval(name))){
    print("ACHIEVE ASSOCIATION FAILED")
    return("NO ASSOCIATION SNPs")
  }else{
    print("expDS")
    expDS <- get(eval(name))
  }

  if(is.null(expDS)){
    return("expDS is NULL")
  }
  
  temp_c<-merge(expDS,pcos_outcome_dat,by.x="SNP",by.y="SNP")
  if(is.null(temp_c)){
    return(NULL)
    }else{
    write.csv(temp_c,paste0("term_PCOS/",id1Nam,"_PCOS.csv"))
    chd_out_dat<-read_outcome_data(
      snps=expDS$SNP,
      filename= paste0("term_PCOS/",id1Nam,"_PCOS.csv"),
      sep=',',
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      eaf_col = "eaf.outcome",
      pval_col = "pval.outcome"
    )
  }
  print('chd_our_dat has been generated')
  
  
  print("harmonise")
  # harmonise data
  dat <- tryCatch(
    harmonise_data(
      exposure_dat = expDS, 
      outcome_dat = chd_out_dat
    )
    ,
    error=function(e){
      cat("Error in harmonise_data. Skipping to the next iteration.\n")
      return(NULL)
      
    }
  )
  if (is.null(dat)) {
    cat("Skipping analysis due to missing data.\n")
    return(NULL)
  }
  
  
  #analysis 
  res<-mr(dat)
  print('res finish')
  
  if(nrow(res)==0){
    print("RES 0 ROW")
    # stop("RES 0 ROW")
    return("RES 0 ROW")
  }
  
  # Tidy up results and save
  res_tidy<-res %>%
    split_outcome() %>% 
    split_exposure() %>%
    generate_odds_ratios()
  print("analysis finish")
  # return(res_tidy)
  res_tidy$'factor'<-id1Cat
  assign(paste0("Res_tidy2_",id1Nam), res_tidy, envir = .GlobalEnv)
  
  print("Heterogeneity and Pleiotropy")
  if(nrow(dat)<=3){
    return("dat rows <= 3 ")
  }
  mr_heterogeneity(dat,method_list=c("mr_ivw_mre","mr_egger_regression"))
  out_sensitivity <- full_join(mr_pleiotropy_test(dat), 
                               mr_heterogeneity(dat)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list <- list(dataset1 = res_tidy, dataset2 = out_sensitivity)
  assign(paste0("Out_sen2_",id1Nam), out_sensitivity, envir = .GlobalEnv)
  
  print("Sensitivity analysis")
  # sensitivity analysis
  sen <- mr_leaveoneout(dat)
  
  
  ###########################################PCOS AS EXPOSURE#########################################
  print("START ANALYSIS PCOS AS EXPOSURE")
  chd_out_dat2 <- tryCatch(
    extract_outcome_data(
      snps = pcos_exp_dat2$SNP,
      outcomes = id1
    ),
    error=function(e){
      print("chd_out_dat ERROR")
    }
  )
  
  if (is.null(chd_out_dat2)) {
    cat("Skipping analysis due to missing chd_out_dat.\n")
    return(NULL)
  }
  
  print("chd_out_dat2 has been genereated")
  
  dat2 <- harmonise_data(
    exposure_dat = pcos_exp_dat2,
    outcome_dat = chd_out_dat2
  )

  if(nrow(dat2)==0){
    return("Do not have enough dat2 for analysis")
  }
  print("Start to anaysis")
  
  res2<-tryCatch(
    mr(dat2),
    error=function(e){
      print("analysis ERROR")
    }
  )
  
  if(nrow(res2)==0){
    print("RES 0 ROW")
    return("RES 0 ROW")
  }
  
  res_tidy2<-res2 %>%
    split_outcome() %>%
    split_exposure() %>%
    generate_odds_ratios()
  res_tidy2$'factor'<-id1Cat
  assign(paste0("Bi_Res_tidy2_",id1Nam), res_tidy2, envir = .GlobalEnv)
  
  
  if(nrow(dat2)<=3){
    return("dat2 rows <= 3 ")
  }
  out_sensitivity2 <- full_join(mr_pleiotropy_test(dat2),
                                mr_heterogeneity(dat2)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list2 <- list(dataset1 = res_tidy2, dataset2 = out_sensitivity2)
  assign(paste0("Bi_Out_sen2_",id1Nam), out_sensitivity2, envir = .GlobalEnv)
  
  print('FINISH BIDIRECTION ANALYSIS')
  
}


#ANALYSIS FOR BIDIRECTION
bi_ana_bat3<-function(id1,id1Nam,id1Cat){
  # id1<-"ebi-a-GCST90013974"
  # id1Nam<-38
  # id1Cat<-"Physical Trait"
  ###########################################PCOS AS OUTCOME#########################################
  # exact the exposure from the IEU DATABASE
  print("START ANALYSIS PCOS AS OUTCOME")
  name <- paste0('exp_dat_',id1Nam)
  assign(eval(name),extract_instruments(eval(id1)),envir = .GlobalEnv)
  
  if(!exists(eval(name))){
    print("ACHIEVE ASSOCIATION FAILED")
    return("NO ASSOCIATION SNPs")
  }else{
    print("expDS")
    expDS <- get(eval(name))
  }
  
  print('chd_our_dat')
  # PCOS as outcomes
  # chd_out_dat <- extract_outcome_data(
  #   snps = expDS$SNP,
  #   outcomes = "finn-b-E4_POCS"
  # )
  chd_out_dat <- tryCatch(
    extract_outcome_data(
      snps = expDS$SNP,
      outcomes = "ebi-a-GCST90044902"
    ),
    error=function(e){
      # print("chd_out_dat ERROR")
      stop("chd_out_dat ERROR")
    }
  )
  
  if (is.null(chd_out_dat)) {
    cat("Skipping analysis due to missing chd_out_dat.\n")
    return(NULL)
  }
  
  
  print("harmonise")
  # harmonise data
  dat <- tryCatch(
    harmonise_data(
      exposure_dat = expDS, 
      outcome_dat = chd_out_dat
    )
    ,
    error=function(e){
      cat("Error in harmonise_data. Skipping to the next iteration.\n")
      return(NULL)
      
    }
  )
  if (is.null(dat)) {
    cat("Skipping analysis due to missing data.\n")
    return(NULL)
  }
  
  
  #analysis 
  res<-mr(dat)
  print('res finish')
  
  if(nrow(res)==0){
    print("RES 0 ROW")
    # stop("RES 0 ROW")
    return("RES 0 ROW")
  }
  
  # Tidy up results and save
  res_tidy<-res %>%
    split_outcome() %>% 
    split_exposure() %>%
    generate_odds_ratios()
  print("analysis finish")
  # return(res_tidy)
  res_tidy$'factor'<-id1Cat
  assign(paste0("Res_tidy3_",id1Nam), res_tidy, envir = .GlobalEnv)
  
  if(nrow(dat)<=3){
    return("dat rows <= 3 ")
  }
  
  print("Heterogeneity and Pleiotropy")
  mr_heterogeneity(dat,method_list=c("mr_ivw_mre","mr_egger_regression"))
  out_sensitivity <- full_join(mr_pleiotropy_test(dat), 
                               mr_heterogeneity(dat)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list <- list(dataset1 = res_tidy, dataset2 = out_sensitivity)
  assign(paste0("Out_sen3_",id1Nam), out_sensitivity, envir = .GlobalEnv)
  
  print("Sensitivity analysis")
  # sensitivity analysis
  sen <- mr_leaveoneout(dat)
  
  
  ###########################################PCOS AS EXPOSURE#########################################
  print("START ANALYSIS PCOS AS EXPOSURE")
  
  chd_out_dat2 <- tryCatch(
    extract_outcome_data(
      snps = pcos_exp_dat3$SNP,
      outcomes = id1
    ),
    error=function(e){
      print("chd_out_dat ERROR")
      # stop("chd_out_dat ERROR")
    }
  )
  
  if (is.null(chd_out_dat2)) {
    cat("Skipping analysis due to missing chd_out_dat.\n")
    return(NULL)
  }
  
  dat2 <- harmonise_data(
    exposure_dat = pcos_exp_dat3,
    outcome_dat = chd_out_dat2
  )
  if (is.null(dat2)) {
    cat("Skipping analysis due to dat2.\n")
    return(NULL)
  }
  
  res2<-mr(dat2)
  res_tidy2<-res2 %>%
    split_outcome() %>%
    split_exposure() %>%
    generate_odds_ratios()
  res_tidy2$'factor'<-id1Cat
  assign(paste0("Bi_Res_tidy3_",id1Nam), res_tidy2, envir = .GlobalEnv)
  
  
  
  if(nrow(dat2)<=3){
    return("dat2 rows <= 3 ")
  }
  out_sensitivity2 <- full_join(mr_pleiotropy_test(dat2),
                                mr_heterogeneity(dat2)) %>%
    dplyr::select(-c("id.outcome", "id.exposure"))
  # res_full<-full_join(res_tidy, out_sensitivity, by=c('id.exposure','method'))
  datasets_list2 <- list(dataset1 = res_tidy2, dataset2 = out_sensitivity2)
  assign(paste0("Bi_Out_sen3_",id1Nam), out_sensitivity2, envir = .GlobalEnv)
  
  print('FINISH BIDIRECTION ANALYSIS')
  
}


loop_append<-function(list){
  combined_temp<-NULL
  for (obj_name in list) {
    if(is.list(get(obj_name))){
      ds<-get(obj_name)
      ds$'resource' <- obj_name
      if(is.null(combined_temp)){
        print("Initial")
        combined_temp<-ds
      }else{
        print(obj_name)
        combined_temp <- rbind(combined_temp, ds)
      }
    }
  }
 return(combined_temp)
}
  




newanalysis<-function(id1,id1Nam,id1type){
  # id1<-"ieu-b-4870"
  # id1Nam<-"SHBG"
  # id1type<-"continuous"
  name <- paste0('exp_dat_',id1Nam)
  if(exists(name)){
    print(paste(name,"has been exists"))
  }else{
    assign(eval(name),extract_instruments(eval(id1)),envir = .GlobalEnv)
  }
  
  if(!exists(eval(name))){
    print("ACHIEVE ASSOCIATION FAILED")
    return("NO ASSOCIATION SNPs")
  }else{
    print("expDS")
    expDS <- get(eval(name))
  }
 
  
  chd_out_dat <-extract_outcome_data(
    snps = expDS$SNP,
    outcomes = "ebi-a-GCST90044902"
  )
  dat<-harmonise_data(
    exposure_dat = expDS, 
    outcome_dat = chd_out_dat
  )
  print("MRMixx")
  
  data = dat %>% mutate(MAF = pmin(eaf.exposure,1-eaf.exposure))

  data_std<-tryCatch(
    with(data, standardize(beta.exposure,beta.outcome,se.exposure,se.outcome,xtype = id1type, ytype = "binary", nx = samplesize.exposure, ny = NULL, MAF = MAF))
    ,
    error=function(e){
      return(NULL)
    }
  )
  if (is.null(data_std)) {
    cat("Skipping analysis due to missing data_std \n")
    return("NULL")
  }

  res1<- MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std, profile = TRUE)
  res1_result<-data.frame(res1[c(1,4,6)])
  colnames(res1_result)<-c('b.Mixx','se.Mixx','pval.Mixx')
  res2<-mr(dat,method="mr_ivw")[c(4,7:9)]
  combined<-cbind(res2, res1_result)
  assign(paste0("Check_",id1Nam), combined, envir = .GlobalEnv)
  
}


multiple_ana<-function(explist,outid){
  id_outcome <- outid
  id_exposure <- explist
  exposure_dat<- mv_extract_exposures(id_exposure)
  oucome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
  if(is.null(oucome_dat)){
    return("Stop for missing data in oucome_dat")
  }
  mvdat<- mv_harmonise_data(exposure_dat, oucome_dat)  
  mvres<- mv_multiple(mvdat) 
  # mvres<-mv_ivw(mvdat)
  return(mvres)
}

