
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(TMB)
library(optparse)
library(gridExtra)

enough_samples = read.table("../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]

source("1_create_ROO/roo_functions.R")
source("2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("2_inference_TMB/mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")
source("../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

opt = list()

re_run <- F

opt$output = NA
opt$model = "fullREDM" 
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

flder <- "../data/roo/"
fles <- list.files(flder)
fles <- fles[grepl('signaturesPCAWG', fles)]
fles <- fles[gsub("_signaturesPCAWG_ROO.RDS", "", fles) %in% enough_samples]
# opt$input = "../data/roo/Uterus-AdenoCA_signaturesPCAWG_ROO.RDS"

for(xx in fles){
  
  try(rm(dataset_amalgamated))
  try(rm(dataset_subset))
  try(rm(results_inference0))
  
  opt$input <- paste0(flder, xx)
  name_plots <- gsub("_ROO.RDS", "", basename(opt$input))
  
  cat(opt$input)
  dataset = load_PCAWG(ct = opt$input, typedata = opt$feature_type, simulation = simulation_bool,
                       path_to_data = NA, read_directly=opt$read_directly)
  
  dataset_roo <- readRDS(opt$input)
  dataset_roo_new <- dataset_roo
  # give_barplot_from_obj(dataset)
  
  if(opt$nonexo_bool | grepl('nonexo', opt$output)){
    ## select only nonexogenous signatures
    nonexogenous = read.table("../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                              comment.char = "#", fill = F)
    dataset <- give_subset_sigs_TMBobj(dataset, sigs_to_remove = nonexogenous$V1)
  }
  
  # give_barplot_from_obj(dataset)
  
  if(opt$nonexo_bool){
    results_inference <- readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnonexo_", name_plots, ".RDS" ))
  }else{
    results_inference <- readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", name_plots, ".RDS" ))
  }
  plot_betas(results_inference)
  # results_inference0 = try(wrapper_run_TMB(object = dataset, model = mod_model_name, use_nlminb=use_nlminb))
  # saveRDS(paste0("../../../../data/pcawg_robjects_cache/fullREDMnoscalingnonexo_convergence/", name_plots, '_', mod_model_name, '.RDS'))
  
  dmin1 <- length(python_like_select_name(results_inference$par.fixed, 'beta'))/2
  
  if(dmin1 > 1){
    ## not correct - cov_par_RE are the values of L, not the off-entries of the covariance matrix
    cov_estimates <- fill_covariance_matrix(arg_d = dmin1,
                                            arg_entries_var = rep(1, dmin1),
                                            arg_entries_cov = python_like_select_name(results_inference$par.fixed, 'cov_par_RE'))
    colnames(cov_estimates) <- rownames(cov_estimates) <- vector_cats_to_logR(colnames(dataset$Y))
    
    # results_inferenc_diagDMDL = try(wrapper_run_TMB(object = dataset, model = "diagRE_DM", use_nlminb=use_nlminb))
    # results_inference$pdHess
    # results_inferenc_diagDMDL$pdHess
    
    pdf(paste0("../results/models_explanatory/fullREDMnoscalingnonexo_convergence/", name_plots, '_estimated_cov_matrix.pdf'))
    print(pheatmap::pheatmap(cov_estimates, main = paste0(name_plots, '\nEstimates of covariances')))
    dev.off()
    
    # pdf(paste0("../results/models_explanatory/fullREDMnoscalingnonexo_convergence/", name_plots, '_all_cov_of_estimates.pdf'))
    # print(pheatmap::pheatmap(results_inference$cov.fixed, main = name_plots))
    # dev.off()
    # results_inference
    
    ## what are those signatures with super high covariances?
    report <- TMB::summary.sdreport(results_inference)
    
    names_covariances <- outer(rownames(cov_estimates), rownames(cov_estimates), Vectorize(function(i,j) {paste0(i,'-', j)}))
    names_covariances_flat <- names_covariances[lower.tri(names_covariances)]
    names_covariances[upper.tri(names_covariances)] ## same but inverted (i.e. sign will be opposite, because they are covariances), but same order
    
    colnames(dataset$Y)[-ncol(dataset$Y)]
    rownames(results_inference$cov.fixed)[grepl('cov', rownames(results_inference$cov.fixed))] <- paste0('cov_', names_covariances_flat) #paste0(colnames(dataset$Y)[-ncol(dataset$Y)],
    #      rownames(results_inference$cov.fixed)[grepl('cov',
    #                                                  rownames(results_inference$cov.fixed))])
    colnames(results_inference$cov.fixed)[grepl('cov', colnames(results_inference$cov.fixed))] <-  paste0('cov_', names_covariances_flat) #paste0(colnames(dataset$Y)[-ncol(dataset$Y)],
    #        colnames(results_inference$cov.fixed)[grepl('cov',
    #                                                    colnames(results_inference$cov.fixed))])
    pdf(paste0("../results/models_explanatory/fullREDMnoscalingnonexo_convergence/", name_plots, '_all_cov_annotated.pdf'), width = 10, height = 10)
    print(pheatmap::pheatmap(results_inference$cov.fixed, main = paste0(name_plots, '\nCovariances of estimates of covariances')))
    dev.off()
    
    
    subcovmat <- results_inference$cov.fixed[grepl('cov', rownames(results_inference$cov.fixed)),
                                             grepl('cov', colnames(results_inference$cov.fixed))]
    pdf(paste0("../results/models_explanatory/fullREDMnoscalingnonexo_convergence/", name_plots, '_all_cov_covRE.pdf'), width = 10, height = 10)
    if(is.null(dim(subcovmat))){
      plot(0,0)
    }else{
      print(pheatmap::pheatmap(subcovmat, main = paste0(name_plots, '\nCovariances of estimates of covariances')))
    }
    dev.off()
    
    ## -------- subset of signatures  
    results_inference$pdHess
    
    if(xx == "Bone-Osteosarc_signaturesPCAWG_ROO.RDS"){
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(c('SBS17a', 'SBS17b')),
                                                                                                      as.list(colnames(dataset$Y)[!grepl('SBS17', colnames(dataset$Y))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))
      
    }else if (xx == "ColoRect-AdenoCA_signaturesPCAWG_ROO.RDS"){
      
      colMeans(normalise_rw(dataset$Y))
      
      dataset_subset <- give_subset_sigs_TMBobj(sig_obj = dataset,  sigs_to_remove = c('SBS17a'))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_subset, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_subset$Y[dataset_subset$x[,2] == 0,],
                                                dataset_subset$Y[dataset_subset$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))      
      
    }else if(xx == "Head-SCC_signaturesPCAWG_ROO.RDS"){
      
      colMeans(normalise_rw(dataset$Y))
      
      amalgamation_vector <- c('SBS17a', 'SBS17b')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) grep(i, colnames(dataset$Y)))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))      
      
      
    }else if(xx == "Kidney-ChRCC_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      amalgamation_vector <- c('SBS17a', 'SBS17b')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) grep(i, colnames(dataset$Y)))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }else if(xx == "Kidney-RCC.papillary_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      amalgamation_vector <- c('SBS40', 'SBS41')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) grep(i, colnames(dataset$Y)))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }else if(xx == "Lymph-BNHL_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      pairs(normalise_rw(dataset$Y))
      
      dataset_subset <- give_subset_sigs_TMBobj(sig_obj = dataset,  sigs_to_remove = c('SBS13'))
      # results_inference0 = try(wrapper_run_TMB(object = dataset_subset, model = mod_model_name, use_nlminb=use_nlminb))
      # results_inference0$pdHess
      
      amalgamation_vector <- c('SBS17a', 'SBS17b')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_subset,  list_groupings = c(list(amalgamation_vector),
                                                                                                             as.list(colnames(dataset_subset$Y)[-sapply(amalgamation_vector, function(i) grep(i, colnames(dataset_subset$Y)))])))
      amalgamation_vector <- c('SBS6', 'SBS34', 'SBS40')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_amalgamated,  list_groupings = c(list(amalgamation_vector),
                                                                                                                  as.list(colnames(dataset_amalgamated$Y)[-sapply(amalgamation_vector, function(i) grep(i, colnames(dataset_amalgamated$Y)))])))
      pairs(dataset_amalgamated$Y)
      
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }else if(xx == "Ovary-AdenoCA_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      
      amalgamation_vector <- c('SBS2', 'SBS40')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset$Y) == i))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
      
    }else if(xx == "Panc-Endocrine_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      # pairs(normalise_rw(dataset$Y))
      amalgamation_vector <- c('SBS26', 'SBS30', 'SBS36', 'SBS39')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset$Y) == i))])))
      pairs(normalise_rw(dataset_amalgamated$Y))
      # amalgamation_vector <- c('SBS2', 'SBS13')
      # dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_amalgamated,  list_groupings = c(list(amalgamation_vector),
      # as.list(colnames(dataset_amalgamated$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset_amalgamated$Y) == i))])))
      pairs(normalise_rw(dataset_amalgamated$Y))
      # dataset_subset <- give_subset_sigs_TMBobj(sig_obj = dataset,  sigs_to_remove = c('SBS2'))
      
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }else if(xx == "Prost-AdenoCA_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      amalgamation_vector <- c('SBS2', 'SBS40')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset$Y) == i))])))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_amalgamated, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
      
    }else if(xx == "Stomach-AdenoCA_signaturesPCAWG_ROO.RDS"){
      colMeans(normalise_rw(dataset$Y))
      pairs(normalise_rw(dataset$Y))
      amalgamation_vector <- c('SBS17a', 'SBS17b')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset,  list_groupings = c(list(amalgamation_vector),
                                                                                                      as.list(colnames(dataset$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset$Y) == i))])))
      amalgamation_vector <- c('SBS2', 'SBS13', 'SBS40', 'SBS41')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_amalgamated,  list_groupings = c(list(amalgamation_vector),
                                                                                                                  as.list(colnames(dataset_amalgamated$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset_amalgamated$Y) == i))])))
      # pairs(normalise_rw(dataset_amalgamated$Y))
      dataset_subset <- give_subset_sigs_TMBobj(sig_obj = dataset_amalgamated,  sigs_to_remove = c('SBS28'))
      
      
      pairs(normalise_rw(dataset_subset$Y))
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = dataset_subset, model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      dataset_roo_new@count_matrices_all = list(dataset_subset$Y[dataset_subset$x[,2] == 0,],
                                                dataset_subset$Y[dataset_subset$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }else if(xx == "Uterus-AdenoCA_signaturesPCAWG_ROO.RDS" ){
      colMeans(normalise_rw(dataset$Y))
      dataset_subset <- give_subset_sigs_TMBobj(sig_obj = dataset,  sigs_to_remove = c('SBS28'))
      pairs(normalise_rw(dataset_subset$Y))
      
      ## interestingly, SBS6 did not seem to be problematic but it was
      amalgamation_vector <- c('SBS2', 'SBS13', 'SBS6')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_subset,  list_groupings = c(list(amalgamation_vector),
                                                                                                             as.list(colnames(dataset_subset$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset_subset$Y) == i))])))
      amalgamation_vector <- c('SBS10a', 'SBS10b')
      dataset_amalgamated <- give_amalgamated_exposures_TMBobj(sig_obj = dataset_amalgamated,  list_groupings = c(list(amalgamation_vector),
                                                                                                                  as.list(colnames(dataset_amalgamated$Y)[-sapply(amalgamation_vector, function(i) which(colnames(dataset_amalgamated$Y) == i))])))
      pairs(normalise_rw(dataset_amalgamated$Y))
      
      if(re_run){
        results_inference0 = try(wrapper_run_TMB(object = sort_columns_TMB(dataset_amalgamated), model = mod_model_name, use_nlminb=use_nlminb))
        results_inference0$pdHess
      }
      
      dataset_roo_new@count_matrices_all = list(dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 0,],
                                                dataset_amalgamated$Y[dataset_amalgamated$x[,2] == 1,])
      dataset_roo_new@count_matrices_active <- dataset_roo_new@count_matrices_all
      dataset_roo_new@sample_names <- rownames(dataset_roo_new@count_matrices_all[[1]])
      saveRDS(dataset_roo_new, gsub("signaturesPCAWG", "signaturesPCAWGSaA", opt$input))   
      
    }
    
    list.files("../data/pcawg_robjects_cache/tmb_results/nlminb/")[grepl('fullREDMnoscalingnonexosubset_', list.files("../data/pcawg_robjects_cache/tmb_results/nlminb/"))]
    
    if(re_run)  saveRDS(results_inference0, paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" ))
    # results_inference00 <- readRDS(paste0("../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", name_plots, ".RDS" ))
    # results_inference00
    # 
    # ## compare betas
    # betas1 <- plot_betas(results_inference)
    # betas2 <- plot_betas(results_inference00)
    # grid.arrange(betas1, betas2, nrow=1)
    
  }else{
    pdf(paste0("../results/models_explanatory/fullREDMnoscalingnonexo_convergence/", name_plots, '_only_two_cats.pdf'), width = 10, height = 10)
    plot(0,0)
    dev.off()
  }
  
}

## check wald test results depending on the order
dataset_roo

results_inferenc_diagDMDL1 = try(wrapper_run_TMB(object = dataset, model = "diagRE_DM", use_nlminb=use_nlminb))
results_inferenc_diagDMDL1

dataset_reorder <- dataset
dataset_reorder2 <- dataset
dataset_reorder$Y <- dataset_reorder$Y[,ncol(dataset_reorder$Y):1]
dataset_reorder2$Y <- cbind(dataset_reorder$Y[,-4],dataset_reorder$Y[,4])

results_inferenc_diagDMDL2 = try(wrapper_run_TMB(object = dataset_reorder, model = "diagRE_DM", use_nlminb=use_nlminb))
results_inferenc_diagDMDL2

results_inferenc_diagDMDL3 = try(wrapper_run_TMB(object = dataset_reorder2, model = "diagRE_DM", use_nlminb=use_nlminb))
results_inferenc_diagDMDL3

wald_TMB_wrapper(results_inferenc_diagDMDL1)
wald_TMB_wrapper(results_inferenc_diagDMDL2)
wald_TMB_wrapper(results_inferenc_diagDMDL3)

plot(sort(softmax(c(select_slope_2(python_like_select_name(results_inferenc_diagDMDL1$par.fixed, 'beta')), 0))),
sort(softmax(c((select_slope_2(python_like_select_name(results_inferenc_diagDMDL2$par.fixed, 'beta'))), 0))))
abline(coef = c(0,1), lty='dashed')

pairs(cbind(sort(softmax(c(select_slope_2(python_like_select_name(results_inferenc_diagDMDL1$par.fixed, 'beta')), 0))),
            sort(softmax(c((select_slope_2(python_like_select_name(results_inferenc_diagDMDL2$par.fixed, 'beta'))), 0))),
            sort(softmax(c((select_slope_2(python_like_select_name(results_inferenc_diagDMDL3$par.fixed, 'beta'))), 0)))))
abline(coef = c(0,1), lty='dashed')

