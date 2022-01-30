## simulate COSMIC signatures, then re-extractm and compare

# rm(list=ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

give_plot_bleeding <- function(names_sigs, abundances=NULL, rel_path='../../../', scales_x_free=F, return_dataframe=F,
                               resort_sig_labels=T, exposures_input=NULL){

  ## rel_path: rel_path to main folder for GlobalDA
  
  library(ggplot2)
  source(paste0(rel_path, "code/1_create_ROO/helper_1_create_ROO.R"))
  
  n <- 200
  nlambda <- rpois(n, 260) ## number of mutations per observation
  # names_sigs <- c('SBS1', 'SBS40', 'SBS3', 'SBS17a', 'SBS2')
  k <- length(names_sigs)
  
  sigs_cosmic <- read.table(paste0(rel_path, "data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                            stringsAsFactors = FALSE, sep = ',', header = TRUE)
  rownames(sigs_cosmic) <- apply(sigs_cosmic, 1, function(i) paste0(substr(i[2], 1, 1), '[', i[1], ']', substr(i[2], 3, 3)))
  sigs_cosmic <- sigs_cosmic[,-c(1:2)]
  
  ## get exposures
  if(is.null(abundances)){
    exposures <- MCMCpack::rdirichlet(n=n, rep(1/k, k))
  }else{
    if(!is.null(exposures_input)){
      exposures <- exposures_input
    }else{
      exposures <- MCMCpack::rdirichlet(n=n, abundances*5)
    }
  }
  
  nlambda_sig <- t(sapply(1:n, function(n_it){
    # if(exposures[n_it,] == 0){ ##20220113
    #   rep(0, k)
    # }else{
      rowSums(rmultinom(n = nlambda[n_it], size=1, prob = exposures[n_it,]) )
    # }
    })) ## number of mutations per observation per signature
  colnames(nlambda_sig) <- names_sigs
  
  all(rowSums(nlambda_sig) == nlambda)
  
  ## get mutations
  ## for each signature, choose at random a mutation it creates
  count_obs <- t(apply(nlambda_sig, 1, function(i) {
    rowSums(sapply(1:k, function(i_sig) rowSums(rmultinom(n = i[i_sig], size=1, prob=sigs_cosmic[,names_sigs[i_sig]] )  )))
    }))
  
  all(rowSums(count_obs) == nlambda)
  min(nlambda)
  
  colnames(count_obs) <- rownames(sigs_cosmic)
  rownames(count_obs) <- paste0('Sample', 1:n)
  
  ## get exposures (recovery)
  count_obs_sample = count_obs[1,]
  sigsQP <- t(apply(count_obs, 1, function(count_obs_sample) sum(count_obs_sample)*QPsig(tumour.ref = rep(names(count_obs_sample), count_obs_sample),
                                                         signatures.ref = as(sigs_cosmic[,names_sigs], 'matrix'))))
  
  df_compare_bleeding <- cbind.data.frame(exposuressigsQP=round(as.vector(sigsQP)),
                                          exposures=as.vector(nlambda_sig),
                                          sig=rep(names_sigs, each=n),
                                          sample=rep(rownames(count_obs), k))#,
                                          #expSBS3=rep(nlambda_sig[,'SBS3'], k))
  if(return_dataframe){
    return(df_compare_bleeding)
  }else{
    if(resort_sig_labels)    df_compare_bleeding$sig <- factor(df_compare_bleeding$sig, levels=paste0('SBS', gtools::mixedsort(gsub('SBS', '', unique(df_compare_bleeding$sig)))))
    .a <- ggplot(df_compare_bleeding, aes(x=exposures, y=exposuressigsQP, col=exposures>exposuressigsQP))+
      geom_abline(slope = 1, intercept = 0, lty='dashed', col='blue')+
      geom_point(alpha=0.7)+theme_bw()+guides(col=F)
    if(scales_x_free){
      .a <- .a+facet_wrap(.~sig, scales = "free_x")
    }else{
      .a <- .a+facet_wrap(.~sig)
    }
    .a
  }
}

