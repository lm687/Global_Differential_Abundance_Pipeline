## Creating dataset to assess the performance of models for parameter recovery

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(TMB)
library(optparse)
library(gridExtra)
library(cowplot)

enough_samples = read.table("../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]

source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")
source("../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

opt = list()

opt$output = NA
opt$model = "fullREDMnoscalingnonexo" 
opt$feature_type = "signaturesPCAWG" 
opt$optimiser = "nlminb" 
opt$simulation_bool = F 
opt$read_directly = T 
opt$use_previous_run_startingvals = T

simulation_bool = opt$simulation_bool
use_nlminb = T

opt$nonexo_bool=T
if(opt$nonexo_bool | grepl('nonexo', opt$output)){
  opt$model <- gsub("nonexo", "", opt$model)
}

cat('Model:', opt$model, '\n')
cat('Feature type:', opt$feature_type, '\n')
cat('Using nlminb:', use_nlminb, '\n')
cat('Simulation boolean:', simulation_bool, '\n')

if(opt$nonexo_bool | grepl('nonexo', opt$output)){
  opt$model <- gsub("nonexo", "", opt$model)
}

# if(opt$model == "fullREM"){
TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
mod_model_name = "fullRE_M"
# }else if(opt$model == "diagREM"){
TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial"))
mod_model_name = "diagRE_M"
# }else if(opt$model == "fullREDM"){
TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
mod_model_name = "fullRE_DM"
# }else if(opt$model == "diagREDM"){
TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
mod_model_name = "diagRE_DM"
# }else if(opt$model =="fullREDMsinglelambda"){
TMB::compile("2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
mod_model_name = "fullREDMsinglelambda"
# use_nlminb=T
# }else if(opt$model =="diagREDMsinglelambda"){
TMB::compile("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
mod_model_name = "diagREDMsinglelambda"
# use_nlminb=T
# }else if(opt$model =="FEDMsinglelambda"){
TMB::compile("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda"))
mod_model_name = "FEDMsinglelambda"
# }else if(opt$model =="fullREDMnoscaling"){
TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling.cpp",  "-std=gnu++17")
dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling"))
mod_model_name = "fullREDMnoscaling"
# }else{
#   stop('Specifiy a valid <model>')
# }

nonexogenous = read.table("../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                          comment.char = "#", fill = F)

flder <- "../data/roo/"
fles <- list.files(flder)
fles <- fles[grepl('signaturesPCAWG', fles)]
fles <- fles[gsub("_signaturesPCAWG_ROO.RDS", "", fles) %in% enough_samples]
# opt$input = "../data/roo/Uterus-AdenoCA_signaturesPCAWG_ROO.RDS" 

xx = fles[1]
xx

pdf("../results/results_TMB/pcawg/betas_and_overdispersion.pdf", height = 9, width = 12)
for(xx in fles){

  opt$input <- paste0(flder, xx)
  name_plots <- gsub("_ROO.RDS", "", basename(opt$input))
  name_plots
  
  
  results_inference_diagDMDL = readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", name_plots, ".RDS" ))
  fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" )
  if(file.exists(fle_fullDMDL)){
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }else{
    fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", name_plots, ".RDS" )
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }
  results_inference_fullDMSL =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", name_plots, ".RDS" ))
  results_inference_fullDMDLonefixedslambda =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", name_plots, ".RDS" ))
  results_inference_fullDMDLonefixedslambda_SaA =     try(readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_",
                                                                         gsub("PCAWG", "PCAWGSaA", name_plots), ".RDS" )))
  if(typeof(results_inference_fullDMDLonefixedslambda_SaA) == "character") results_inference_fullDMDLonefixedslambda_SaA = results_inference_fullDMDLonefixedslambda
  results_inference_diagDMDL$pdHess
  results_inference_fullDMDLnoscaling$pdHess
  results_inference_fullDMSL$pdHess
  results_inference_fullDMDLonefixedslambda$pdHess
  
  colnames_nonexo_notsorted_SP <- colnames(give_subset_sigs_TMBobj(load_PCAWG(ct = gsub("_signaturesPCAWG", "", name_plots),
                                                                              typedata = "signaturesPCAWG",
                                                                              path_to_data = "../data/"), nonexogenous$V1)$Y)
  logR_nonexo_notsorted_SP <- vector_cats_to_logR(colnames_nonexo_notsorted_SP)
  
  # colnames_nonexo_notsorted_SP_subset_and_amalgamation <- colnames(dataset_amalgamated$Y)
  # logR_nonexo_notsorted_SP_subset_and_amalgamation <- vector_cats_to_logR(colnames_nonexo_notsorted_SP_subset_and_amalgamation)
  
  # length(python_like_select_name(results_inference_fullDMDLnoscaling$par.fixed, 'beta'))/2
  # length(logR_nonexo_notsorted_SP_subset)
  # results_inference0
  
  ## checking the overdispersion parameters
  
  
  
  # pltA <- cowplot::as_grob(grid.arrange(plot_betas(results_inference_diagDMDL, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
  #                                  title = paste0('diagDMDL\npval=', wald_TMB_wrapper(results_inference_diagDMDL))),
  #              plot_betas(results_inference_fullDMSL, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
  #                         title=paste0('fullDMSL\npval= ', wald_TMB_wrapper(results_inference_fullDMDLnoscaling))),
  #              plot_betas(results_inference_fullDMDLonefixedslambda, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
  #                         title=paste0('fullDMDLonefixedslambda\npval= ', wald_TMB_wrapper(results_inference_fullDMDLonefixedslambda))),
  #              plot_betas(results_inference_fullDMDLnoscaling, plot =F, return_plot=T,
  #                         title=paste0('fullDMDLnoscaling\npval=', wald_TMB_wrapper(results_inference_fullDMDLnoscaling))),
  #              nrow=4))
  pltA <- cowplot::plot_grid(plot_betas(results_inference_diagDMDL, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
                                        title = paste0('diagDMDL\npval=', wald_TMB_wrapper(results_inference_diagDMDL))),
                             plot_betas(results_inference_fullDMSL, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
                                        title=paste0('fullDMSL\npval= ', wald_TMB_wrapper(results_inference_fullDMDLnoscaling))),
                             plot_betas(results_inference_fullDMDLonefixedslambda, names_cats = logR_nonexo_notsorted_SP, plot =F, return_plot=T,
                                        title=paste0('fullDMDLonefixedslambda\npval= ', wald_TMB_wrapper(results_inference_fullDMDLonefixedslambda))),
                             plot_betas(results_inference_fullDMDLnoscaling, plot =F, return_plot=T,
                                        title=paste0('fullDMDLnoscaling SaA\npval=', wald_TMB_wrapper(results_inference_fullDMDLnoscaling))),
                             plot_betas(results_inference_fullDMDLonefixedslambda_SaA, plot =F, return_plot=T,
                                        title=paste0('fullDMDLonefixedslambda SaA\npval=', wald_TMB_wrapper(results_inference_fullDMDLonefixedslambda_SaA))),
                             nrow=2)
  
  add_fixed_lambda <- function(df){
    df['name'] <- 'Lambda 2'
    df <- rbind.data.frame(df, data.frame(Estimate=0, `Std..Error`=0, name='Lambda 1'))
    df
  }
  overdispersion_plots <- ggplot(melt(list(diagDMDL=plot_lambdas(results_inference_diagDMDL, return_df = T, pl=F),
                                           fullDMSL=plot_lambdas(results_inference_fullDMSL, return_df = T, pl=F),
                                           fullDMDLonefixedslambda=add_fixed_lambda(plot_lambdas(results_inference_fullDMDLonefixedslambda, return_df = T, pl=F)),
                                           fullDMDLonefixedslambdaSaA=add_fixed_lambda(plot_lambdas(results_inference_fullDMDLonefixedslambda_SaA, return_df = T, pl=F)),
                                           fullDMDLnoscaling=plot_lambdas(results_inference_fullDMDLnoscaling, return_df = T, pl=F)),
              id.vars=c('Estimate', 'Std..Error', 'name')),
         aes(x=name, y=`Estimate`))+
    geom_point()+
    geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw()+
    facet_wrap(.~L1)+ggtitle('Overdispersion lambdas')
  
  logSDRE_plots <- ggplot(melt(list(diagDMDL=plot_estimates_TMB(results_inference_diagDMDL, return_df = T, pl=F, parameter_name="logs_sd_RE", verb=F),
                                    fullDMSL=plot_estimates_TMB(results_inference_fullDMSL, return_df = T, pl=F, parameter_name="logs_sd_RE", verb=F),
                                    fullDMDLonefixedslambda=plot_estimates_TMB(results_inference_fullDMDLonefixedslambda, return_df = T, pl=F,
                                                                               parameter_name="logs_sd_RE", verb=F),
                                    fullDMDLonefixedslambdaSaA=plot_estimates_TMB(results_inference_fullDMDLonefixedslambda_SaA, return_df = T,
                                                                               pl=F, parameter_name="logs_sd_RE", verb=F),
                                    fullDMDLnoscaling=plot_estimates_TMB(results_inference_fullDMDLnoscaling,
                                                                         return_df = T, pl=F, parameter_name="logs_sd_RE", verb=F)),
                                      id.vars=c('Estimate', 'Std..Error', 'name')),
                                 aes(x=name, y=`Estimate`))+
    geom_point()+
    geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw()+
    facet_wrap(.~L1)+ggtitle('logs_sd_RE')
  

  # (grid.arrange(pltA,
  #              grid.arrange(overdispersion_plots, logSDRE_plots, nrow=2), ncol=2, top=name_plots))
  print(cowplot::plot_grid(NULL,
                           cowplot::plot_grid(pltA,
                                              cowplot::plot_grid(overdispersion_plots, logSDRE_plots, nrow=1),
                                              ncol=1, rel_heights = c(5,3)),
                           vjust = 1, labels = name_plots, nrow=2, rel_heights = c(1,12)))
  
}

dev.off()


## Check fitted values in each model, and compare it to observed values

name_plots <- "Lymph-CLL_signaturesPCAWG" #gsub("_ROO.RDS", "", basename(paste0(flder[1], xx)))



remove_lambda_pars <- function(i)i[!grepl("log_lambda", names(i))]
remove_sd_pars <- function(i)i[!grepl("logs_sd_RE", names(i))]
remove_lambda_and_cov_pars <- function(i)i[! (grepl("log_lambda", names(i)) | grepl("cov_par", names(i)) )]

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pdf("../results/results_TMB/pcawg/correlation_parameters.pdf", height = 7)
for(xx in fles){
  
  opt$input <- paste0(flder, xx)
  name_plots <- gsub("_ROO.RDS", "", basename(opt$input))
  
  results_inference_diagDMDL = readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", name_plots, ".RDS" ))
  fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" )
  if(file.exists(fle_fullDMDL)){
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }else{
    fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", name_plots, ".RDS" )
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }
  results_inference_fullDMSL =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", name_plots, ".RDS" ))
  results_inference_fullDMDLonefixedslambda =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", name_plots, ".RDS" ))
  
  results_inference_diagDMDL$pdHess
  results_inference_fullDMDLnoscaling$pdHess
  results_inference_fullDMSL$pdHess
  results_inference_fullDMDLonefixedslambda$pdHess

  
  pairs(cbind.data.frame(fullDMDLonefixedslambda=remove_lambda_and_cov_pars(results_inference_fullDMDLonefixedslambda$par.fixed),
                         fullDMSL=remove_lambda_and_cov_pars(results_inference_fullDMSL$par.fixed),
                         diagDMDL=remove_lambda_and_cov_pars(results_inference_diagDMDL$par.fixed)),
        lower.panel = panel.cor, upper.panel = my_line, main=name_plots)
  
}
dev.off()



pdf("../results/results_TMB/pcawg/randomintercepts.pdf", height = 7)
for(xx in fles){
  fullDMDL_is_with_subset <- NA
  
  opt$input <- paste0(flder, xx)
  name_plots <- gsub("_ROO.RDS", "", basename(opt$input))
  
  results_inference_diagDMDL = readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", name_plots, ".RDS" ))
  fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" )
  if(file.exists(fle_fullDMDL)){
    fullDMDL_is_with_subset <- T
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }else{
    fullDMDL_is_with_subset <- F
    fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", name_plots, ".RDS" )
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }
  results_inference_fullDMSL =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", name_plots, ".RDS" ))
  results_inference_fullDMDLonefixedslambda =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", name_plots, ".RDS" ))
  
  results_inference_diagDMDL$pdHess
  results_inference_fullDMDLnoscaling$pdHess
  results_inference_fullDMSL$pdHess
  results_inference_fullDMDLonefixedslambda$pdHess
  
  if(fullDMDL_is_with_subset){
    pairs(cbind.data.frame(diagDMDL=results_inference_diagDMDL$par.random,
                           fullDMDLnoscaling=rep(0, length(results_inference_diagDMDL$par.random)),
                           fullDMSL=results_inference_fullDMSL$par.random,
                           fullDMDLonefixedslambda=results_inference_fullDMDLonefixedslambda$par.random),
          main=name_plots, pch=19, col=scales::alpha('black', 0.2))
  }else{
    pairs(cbind.data.frame(diagDMDL=results_inference_diagDMDL$par.random,
                           fullDMDLnoscaling=results_inference_fullDMDLnoscaling$par.random,
                           fullDMSL=results_inference_fullDMSL$par.random,
                           fullDMDLonefixedslambda=results_inference_fullDMDLonefixedslambda$par.random),
          main=name_plots, pch=19, col=scales::alpha('black', 0.2))
  }
  
}
dev.off()


pdf("../results/results_TMB/pcawg/fitted_logR.pdf", height = 7)
for(xx in fles){
  fullDMDL_is_with_subset <- NA
  
  opt$input <- paste0(flder, xx)
  name_plots <- gsub("_ROO.RDS", "", basename(opt$input))
  
  results_inference_diagDMDL = readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", name_plots, ".RDS" ))
  fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" )
  if(file.exists(fle_fullDMDL)){
    fullDMDL_is_with_subset <- T
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }else{
    fullDMDL_is_with_subset <- F
    fle_fullDMDL <- paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", name_plots, ".RDS" )
    results_inference_fullDMDLnoscaling =     readRDS(fle_fullDMDL)
  }
  results_inference_fullDMSL =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", name_plots, ".RDS" ))
  results_inference_fullDMDLonefixedslambda =     readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", name_plots, ".RDS" ))
  
  results_inference_diagDMDL$pdHess
  results_inference_fullDMDLnoscaling$pdHess
  results_inference_fullDMSL$pdHess
  results_inference_fullDMDLonefixedslambda$pdHess

  models_it <- list(results_inference_diagDMDL, results_inference_fullDMDLnoscaling, results_inference_fullDMSL,
                    results_inference_fullDMDLonefixedslambda)
  logRmats <- lapply(models_it, function(i){  
    dmin1 <- length(python_like_select_name(i$par.fixed, 'beta'))/2
    re_mat <- re_vector_to_matrix(i$par.random, dmin1)
    ntimes2 = nrow(re_mat)*2
    z_matrix <- (give_z_matrix(ntimes2)/2)
    x_matrix <- cbind(1, c(rep(0, ntimes2/2), rep(1, ntimes2/2)))
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(i$par.fixed, 'beta'), nrow=2)
    return(logRmat)
  })
  names(logRmats) <- c('diagDMDL', 'fullDMDLnoscaling', 'fullDMSL',
                    'fullDMDLonefixedslambda')
  
  if(fullDMDL_is_with_subset){
    pairs(do.call('cbind', lapply(logRmats[-2], as.vector)),
          main=name_plots, pch=19, col=scales::alpha(factor(1+rep(c(0,1), each=nrow(logRmats[[1]])/2)), 0.2))
  }else{
    pairs(do.call('cbind', lapply(logRmats, as.vector)),
          main=name_plots, pch=19,  col=scales::alpha(factor(1+rep(c(0,1), each=nrow(logRmats[[1]])/2)), 0.2))
  }
  
}
dev.off()

