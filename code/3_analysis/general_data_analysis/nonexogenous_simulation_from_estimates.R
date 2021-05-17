rm(list = ls())

library(ggrepel)
library(Ternary)
library(MCMCpack)
library(dplyr)
library(gridExtra)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../../2_inference_TMB/helper_TMB.R")
source("../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

enough_samples = readLines("~/Desktop/CT_sufficient_samples.txt")
ct <- enough_samples[5]
ct
nonexogenous = read.table("../../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t", comment.char = "#", fill = F)

# rm(list = c('fullRE_DM', 'fullRE_DM_ns1', 'fullRE_DM_nonexo', 'diagRE_DM', 'sparseRE_DM'))
# fullRE_M <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_M_", ct, "_signatures.RDS"))
# fullRE_M
# fullRE_DM <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_signatures.RDS"))
# fullRE_DM
# fullRE_DM_ns1 <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_nucleotidesubstitution1.RDS"))
# fullRE_DM_ns1
# fullRE_DM_nonexo_old <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DM_", ct, "_signatures.RDS"))
# fullRE_DM_nonexo_old
# fullRE_DM_nonexo <- readRDS(paste0("../../data/../pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DMSL_", ct, "_signatures.RDS"))
# fullRE_DM_nonexo
# sparseRE_DM_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_nonexo_DMSL_", ct, "_signatures.RDS"))
# sparseRE_DM_nonexo
# diagRE_DM <- readRDS(paste0("../../data/../pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS"))
# diagRE_DM
# sparseRE_DM <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS"))
# sparseRE_DM


# wald_TMB_wrapper(fullRE_M, verbatim = F)
# wald_TMB_wrapper(fullRE_DM_ns1, verbatim = F)
# wald_TMB_wrapper(diagRE_DM)
# wald_TMB_wrapper(sparseRE_DM)
# wald_TMB_wrapper(fullRE_DM_nonexo)


diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
})

fullRE_M_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_M_", ct, "_signatures.RDS")))
})



fullRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_M_", ct, "_signatures.RDS")))
})

diagRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/diagRE_M_", ct, "_signatures.RDS")))
})

sparseRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS")))
})

fullRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DMSL_", ct, "_signatures.RDS")))
})

diagRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_nonexo_DMSL_", ct, "_signatures.RDS")))
})

fullRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_signatures.RDS")))
})

diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
})

sparseRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_nonexo_DMSL_", ct, "_signatures.RDS")))
})

fullRE_DMSL_nonexo[sapply(sparseRE_DMSL_nonexo, typeof) == "character"]

pvals_fullRE_M <- sapply(fullRE_M, function(i) try(wald_TMB_wrapper(i)))
pvals_fullRE_M <- p.adjust(pvals_fullRE_M)
pvals_diagRE_DM <- sapply(diagRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_diagRE_DM <- p.adjust(pvals_diagRE_DM)
pvals_DM <- sapply(sparseRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_DM <- p.adjust(pvals_DM)
pvals_DMnonexo <- sapply(sparseRE_DMSL_nonexo, function(i) try(wald_TMB_wrapper(i)))
pvals_DMnonexo <- p.adjust(pvals_DMnonexo)

sort(pvals_DM)
sort(pvals_DMnonexo)

num_samples <- sapply(enough_samples, function(ct){
    .xx <- try(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1)))
    try(nrow(.xx$Y)/2)
})

effect_size <- sapply(enough_samples, function(ct){
  .xx <- try(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                     sigs_to_remove = unique(nonexogenous$V1)))
  try(sum((normalise_rw(.xx$Y[1:(nrow(.xx$Y)/2),]) - normalise_rw(.xx$Y[(1+nrow(.xx$Y)/2):(nrow(.xx$Y)),]))**2)/(nrow(.xx$Y)/2))
})

length(python_like_select_name(fullRE_DM_nonexo$par.fixed, grep_substring = "beta"))
length(python_like_select_name(fullRE_DM_nonexo_old$par.fixed, grep_substring = "beta"))
length(python_like_select_name(sparseRE_DM_nonexo$par.fixed, grep_substring = "beta"))


ggplot(cbind.data.frame(pvals_DM=pvals_DM, pvals_DMnonexo=pvals_DMnonexo,
                        num_samples=as.numeric(num_samples),
                        ct=enough_samples),
       aes(x=log(pvals_DM), y=log(pvals_DMnonexo), col=num_samples,
           label=ct))+geom_point()+
  geom_hline(yintercept = log(0.05))+geom_vline(xintercept = log(0.05))+
  geom_label_repel()

ggplot(cbind.data.frame(pvals_DM=pvals_DM, pvals_DMnonexo=pvals_DMnonexo,
                        num_samples=as.numeric(num_samples),
                        ct=enough_samples),
       aes(x=(pvals_DM), y=(pvals_DMnonexo), col=num_samples,
           label=ct))+geom_point()+
  geom_hline(yintercept = (0.05))+geom_vline(xintercept = (0.05))+
  geom_label_repel()

plot(pvals_fullRE_M, pvals_DM)

ggplot(cbind.data.frame(effect_size=as.numeric(effect_size), minlogpval=-log(as.numeric(pvals_DM)),
                        label=enough_samples), aes(x=effect_size, y=minlogpval, label=label))+
  geom_point()+geom_label_repel(alpha=0.2)


ct <- "Breast-AdenoCA"

## Only for some CT is the cov matrix positive semi-definite!

ct <- enough_samples[5]
ct <- "Prost-AdenoCA"
subset_sigs_sparse_cov_idx_nonexo <- read.table("../../../current/subset_sigs_sparse_cov_idx_nonexo.txt", stringsAsFactors = F, fill = T)
sparseRE_DMSL_nonexo$`Breast-AdenoCA`

give_sim_from_estimates <- function(ct, typedata = "signatures", sigs_to_remove="", model="sparseRE_DM",
                                    bool_nonexo=TRUE, bool_give_PCA, sig_of_interest='SBS8'){
  
  if(model == "fullRE_M"){
    if(!bool_nonexo)    list_estimates <- fullRE_M
    if(bool_nonexo)    list_estimates <- fullRE_M_nonexo
  }else if(model == "fullRE_DM"){
    if(!bool_nonexo)    list_estimates <- fullRE_DMSL
      if(bool_nonexo)    list_estimates <- fullRE_DMSL_nonexo
    }else if(model == "sparseRE_DM"){
    if(bool_nonexo)    list_estimates <- sparseRE_DMSL_nonexo
  }
  
  obj_data <- sort_columns_TMB(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = typedata, path_to_data = "../../../data/"),
                                         sigs_to_remove = sigs_to_remove))
  dmin1 <- ncol(obj_data$Y)-1
  cov_vec = rep(0, (dmin1**2-dmin1)/2)
  
  if(model %in% c("fullRE_M", "fullRE_DM")){
    cov_vec = python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_par_RE')
    ###**** I AM NOT SURE ABOUT THIS BIT BELOW! ARE THEY SD OR VAR???*****###
    ### implementing them as though they were sd ###
    var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))**2
    var_vec_v2 = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))
  }else if(model == "sparseRE_DM"){
    cov_vec[as.numeric(strsplit(subset_sigs_sparse_cov_idx_nonexo[subset_sigs_sparse_cov_idx_nonexo$V1 == ct,"V2"], ',')[[1]])] = python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_RE_part')
    var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))
  }
  
  cov_mat <- fill_covariance_matrix(arg_d = dmin1,
                                    arg_entries_var = var_vec,
                                    arg_entries_cov = cov_vec)
  cov_mat_v2 <- fill_covariance_matrix(arg_d = dmin1,
                                    arg_entries_var = var_vec_v2,
                                    arg_entries_cov = cov_vec)
  # cov_matb <- fill_covariance_matrix(arg_d = dmin1,
  #                                   arg_entries_var = var_vec**2,
  #                                   arg_entries_cov = cov_vec)

  beta_mat = matrix(python_like_select_name(list_estimates[[ct]]$par.fixed, 'beta'), nrow=2)
  
  n_sim = 1000
  x_sim = cbind(1, rep(c(0,1), n_sim))
  u_sim = mvtnorm::rmvnorm(n = n_sim, mean = rep(0,dmin1), sigma = cov_mat)
  
  theta = x_sim %*% beta_mat + (give_z_matrix(n_sim*2)) %*% u_sim
  
  if(model %in% c('sparseRE_DM', 'fullRE_DM')){
    alpha = softmax(cbind(theta, 0))*exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'log_lambda'))
  }else if(model %in% c("fullRE_M")){
    alpha = softmax(cbind(theta, 0))
  }else{
    stop('Check softmax step')
  }

  if(model %in% c('sparseRE_DM', 'fullRE_DM')){
    probs = t(apply(alpha, 1, MCMCpack::rdirichlet, n=1))
  }else if(model %in% c("fullRE_M")){
    probs = alpha
  }

  probs_obs = normalise_rw(obj_data$Y)
  all_probs = rbind(probs_obs, probs)
  
  if(bool_give_PCA){
    pca <- prcomp(all_probs)
    df_pca <- cbind.data.frame(pca=pca$x[,1:2], col=c(rep('Observed', nrow(probs_obs)), rep('Simulated',nrow(probs))),
                     sig_of_interest=all_probs[,sig_of_interest],
                     group=c('early','late'))
    return(list(df_pca, ggplot(df_pca, aes(x=pca.PC1, y=pca.PC2, col=sig_of_interest))+
                  geom_point(alpha=0.7)+facet_wrap(.~interaction(col,group))))
  }else{
    return(all_probs)
  }
  
}

ct <- enough_samples[10]

pca_sim <- give_sim_from_estimates(ct, "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T)

give_sim_from_estimates(enough_samples[2], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T)
give_sim_from_estimates(enough_samples[3], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS2')
give_sim_from_estimates(enough_samples[4], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')
give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')
# give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), model = "fullRE_M", bool_give_PCA = T, sig_of_interest = 'SBS8')

ggsave(paste0("../../results/results_TMB/pcawg/simulation_from_estimates/", ct, ".pdf"))
# counts = apply(alpha, 1, function(i) 

give_sim_from_estimates(enough_samples[5], "signatures", sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T, sig_of_interest = 'SBS8')

CNS_Medullo <- give_subset_sigs_TMBobj(load_PCAWG(enough_samples[5], "signatures", path_to_data = "../../../data/"),
                                sigs_to_remove = unique(nonexogenous$V1))

give_barplot_from_obj(CNS_Medullo, T)

# mvtnorm::rmvnorm(n = 100, mean = rep(0,dmin1), sigma = cov_matb)



amalgamation1 <- as(normalise_rw(all_probs[,10:12]), 'matrix')
rownames(amalgamation1)[rownames(amalgamation1) == ""] = paste0('Sim', 1:n_sim)
bad_idx <- as.numeric(which( (is.na(rowSums(amalgamation1))) | !(rowSums(amalgamation1) == 1)))
amalgamation1 <- amalgamation1[-bad_idx,]
df_pca_amalgamation1 <- df_pca[-bad_idx,]



aa <- sapply(c('Observed', 'Simulated'), function(arg_group){
  am <- amalgamation1[df_pca$col[-bad_idx] == arg_group,]
  dens <- TernaryDensity(am, resolution = 10L)
  plot.new()
  TernaryPlot(atip = colnames(am)[1], btip = colnames(am)[2], ctip = colnames(am)[3],
              grid.lines = 0, grid.col = NULL)
  ColourTernary(dens, direction = 1)
  # ColourTernary(amalgamation, spectrum=rainbow(256))
  TernaryPoints(am, col = factor(df_pca$col[-bad_idx]), pch = '.', cex=5)
  TernaryDensityContour(am, resolution = 30L)
})  

aa[1]
aa[2]

eigen(cov_mat)$values
eigen(cov_matb)$values

#obj_nonexo$z %*% 

table(sapply(fullRE_M_nonexo, function(i) try(i$pdHess)))
table(sapply(fullRE_M, function(i) try(i$pdHess)))
table(sapply(diagRE_M, function(i) try(i$pdHess)))
table(sapply(diagRE_DMSL, function(i) try(i$pdHess))) ## very good
table(sapply(fullRE_DMSL_nonexo, function(i) try(i$pdHess))) ## half-half
table(sapply(diagRE_DMSL_nonexo, function(i) try(i$pdHess))) ## very good

fullRE_DM_Prost_AdenoCA <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", "Prost-AdenoCA", "_signatures.RDS"))
fullRE_DM_Prost_AdenoCA$pdHess
fullRE_DM_Prost_AdenoCA_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DM_Prost-AdenoCA_signatures.RDS"))
fullRE_DM_Prost_AdenoCA_nonexo$pdHess
fullRE_M_Prost_AdenoCA_nonexo <- readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_M_Prost-AdenoCA_signatures.RDS"))
fullRE_M_Prost_AdenoCA_nonexo$pdHess

give_sim_from_estimates(ct = "Prost-AdenoCA", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS8', model = "fullRE_M")
ct = "Prost-AdenoCA"
typedata = "signatures"
bool_nonexo = T
sigs_to_remove=unique(nonexogenous$V1)
bool_give_PCA = T
sig_of_interest = 'SBS8'
model = "fullRE_M"

fullRE_DM_Prost_AdenoCA

give_sim_from_estimates(ct = enough_samples[2], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")

give_sim_from_estimates(ct = enough_samples[9], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")


give_sim_from_estimates(ct = enough_samples[10], typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")

give_sim_from_estimates(ct = "Lung-SCC", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_DM")

give_sim_from_estimates(ct = "Lung-SCC", typedata = "signatures", bool_nonexo = T,
                        sigs_to_remove=unique(nonexogenous$V1), bool_give_PCA = T,
                        sig_of_interest = 'SBS2', model = "fullRE_M")


fullRE_DMSL_nonexo$`Lymph-CLL`$pdHess
fullRE_M_nonexo$`Lymph-CLL`$pdHess

Lymph_CLL_obj <- give_subset_sigs_TMBobj(load_PCAWG(ct = "Lymph-CLL", typedata = "signatures", path_to_data = "../../../data/"),
                        sigs_to_remove = unique(nonexogenous$V1))
give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                            data_object = Lymph_CLL_obj, print_plot = T, nreps = 20, model = "DMSL", integer_overdispersion_param = 1)

good_plots = give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                                         data_object = Lymph_CLL_obj, print_plot = F, nreps = 20, model = "DMSL", integer_overdispersion_param = 1)
good_plots_unsorted = lapply(list(good_plots), function(i) lapply(i, function(j) cbind.data.frame(sorted_value=as.vector(j), rank_number=1:length(j)) ))
give_interval_plots_2(df_rank = good_plots_unsorted[[1]], data_object = Lymph_CLL_obj,
                      loglog = F, title = 'Lymph_CLL')
good_plots_M = give_ranked_plot_simulation(tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`,
                                         data_object = Lymph_CLL_obj, print_plot = F, nreps = 20, model = "M")
good_plots_unsorted_M = lapply(list(good_plots_M), function(i) lapply(i, function(j) cbind.data.frame(sorted_value=as.vector(j), rank_number=1:length(j)) ))
give_interval_plots_2(df_rank = good_plots_unsorted_M[[1]], data_object = Lymph_CLL_obj,
                      loglog = F, title = 'Lymph_CLL')

pdf("../../../results/results_TMB/pcawg/assessing_models/ranked_plots/Lymph_CLL_DMSL_vs_M.pdf",
    height = 3.5, width = 5)
grid.arrange(give_interval_plots_2(df_rank = good_plots_unsorted_M[[1]], data_object = Lymph_CLL_obj,
                                   loglog = F, title = 'Lymph CLL (M)')+guides(col=guide_legend(title="In confidence interval")),
             give_interval_plots_2(df_rank = good_plots_unsorted[[1]], data_object = Lymph_CLL_obj,
                                   loglog = F, title = 'Lymph CLL (DM)')+guides(col=guide_legend(title="In confidence interval")),
             nrow=1
)
dev.off()

# tmb_fit_object = fullRE_DMSL_nonexo$`Lymph-CLL`
# data_object = Lymph_CLL_obj
# print_plot = T
# nreps = 20
# model = "DMSL"
# integer_overdispersion_param = 1
# 
# unction(tmb_fit_object, data_object, print_plot = T, nreps = 1, model, integer_overdispersion_param){
#   
#   if(model == 'M'){
#     ## theta is always going to be the same. Only replicate the draws from the multinomial
#     
#     sim_theta = simulate_from_M_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                     x_matrix = data_object$x, z_matrix = data_object$z)
#     
#     sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
#     
#     if(nreps>1){
#       sim_counts = replicate(nreps, t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) ), simplify = F)
#     }
#   }else if(model %in% c('DM', 'DMSL')){
#     give_sim_data = function(){
#       if(model == 'DM'){
#         sim_theta = simulate_from_DM_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                          x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
#       }else if(model == 'DMSL'){
#         sim_theta = simulate_from_DMSL_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
#                                            x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
#       }
#       sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
#       return(sim_counts)
#     }    
#     if(nreps>1){
#       sim_counts = replicate(nreps, give_sim_data(), simplify = F)
#     }else{
#       sim_counts = give_sim_data()
#     }
#   }else{
#     stop('Specify a correct model')
#   }
#   
#   stopifnot(all(dim(data_object$Y) == dim(sim_counts)))
#   if(print_plot)  plot(sort(data_object$Y), sort(sim_counts))
#   
#   return(sim_counts)
#   
# }
# 
# 
# tmb_fit_object = tmb_fit_object
# full_RE = T
# x_matrix = data_object$x
# z_matrix = data_object$z
# integer_overdispersion_param=integer_overdispersion_param
# 
# function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
#   dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
#   overdispersion_lambda = rep(integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")), nrow(x_matrix))
#   if(full_RE){
#     re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
#     ntimes2 = nrow(z_matrix)
#     logRmat = z_matrix %*% re_mat + 
#       x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
#     sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
#   }else{
#     sim_thetas = softmax(cbind(sapply(1:dmin1,
#                                       function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
#                                  give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
#     sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
#   }
#   return(sim_thetas)
# }