## similar to posterior_predictive_checks.R, but for the maximum likelihood data

#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggplot2)
require(R.utils)
require(dplyr)
library(parallel)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(MCMCpack) ## sample from Dirichlet
source("../2_inference_TMB/helper_TMB.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("../3_analysis/helper/helper_simulate_from_model.R")
# set.seed(1234)
folder_roo = "../../data/roo/"
folder_robjs = "../../data/robjects_cache/tmb_results/"

coefficient_overdispersion = 1000 ## value used in TMB for better convergence

warning('Be careful to use the latest <coefficient_overdispersion> value!')
#-------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------#
count_objects = sapply(list.files(folder_roo, full.names = TRUE), readRDS)
names(count_objects) = gsub("_ROO.RDS", "", list.files(folder_roo))
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                                   strsplit, split = "_")))
colnames(samples_files) = c('CT', 'type')
table(samples_files[,1], samples_files[,2])
samples_files2 = samples_files %>% filter(type != "nucleotidesubstitution3")
rownames(samples_files2) = rownames(samples_files)[samples_files$type != "nucleotidesubstitution3"]
#-------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------#
results_TMB_M = lapply( python_like_select(list.files(folder_robjs), "^M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_M) = sapply(python_like_select(list.files(folder_robjs), "^M_"), clean_name)
results_TMB_DM = lapply( python_like_select(list.files(folder_robjs), "^DM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_DM) = sapply(python_like_select(list.files(folder_robjs), "^DM_"), clean_name)
results_TMB_DM_dep = lapply( python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"),
                             function(i) readRDS(paste0("../../data/robjects_cache/tmb_results_dep/", i)))
names(results_TMB_DM_dep) = sapply(python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"), clean_name)
results_TMB_LNM = lapply( python_like_select(list.files(folder_robjs), "^LNM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_LNM) = sapply(python_like_select(list.files(folder_robjs), "^LNM_"), clean_name)
results_TMB_fullRE_M = lapply( python_like_select(list.files(folder_robjs), "^fullRE_M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_M) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_M_"), clean_name_fullRE)
full_RE_DM = python_like_select(list.files(folder_robjs), "^fullRE_DM_"); full_RE_DM = full_RE_DM[-grep("_altpar_", full_RE_DM)]
results_TMB_fullRE_DM = lapply( full_RE_DM, function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DM) = sapply(full_RE_DM, clean_name_fullRE)
# results_TMB_fullRE_DM = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DM_altpar_"), function(i) readRDS(paste0(folder_robjs, i)))
# names(results_TMB_fullRE_DM) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DM_altpar_"), clean_name_fullRE_2)
results_TMB_diagRE_M = lapply( python_like_select(list.files(folder_robjs), "^diagRE_DM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_diagRE_M) = sapply(python_like_select(list.files(folder_robjs), "^diagRE_DM_"), clean_name_fullRE)
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
results_TMB_M = results_TMB_M[match(gsub("_", "", names(count_objects)), names(results_TMB_M))]
results_TMB_DM = results_TMB_DM[match(gsub("_", "", names(count_objects)), names(results_TMB_DM))]
results_TMB_LNM = results_TMB_LNM[match(gsub("_", "", names(count_objects)), names(results_TMB_LNM))]
results_TMB_fullRE_M = results_TMB_fullRE_M[match(gsub("_", "", names(count_objects)), names(results_TMB_fullRE_M))]
results_TMB_fullRE_DM = results_TMB_fullRE_DM[match(gsub("_", "", names(count_objects)), names(results_TMB_fullRE_DM))]
results_TMB_diagRE_DM = results_TMB_diagRE_M[match(gsub("_", "", names(count_objects)), names(results_TMB_diagRE_M))]
sapply(list(results_TMB_M, results_TMB_DM, results_TMB_LNM, results_TMB_fullRE_M,
            results_TMB_fullRE_DM, count_objects), length)
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
# M

give_plot_dataset = function(idx_dataset, resultsTMB_list, full_RE=FALSE){
  # if(is.null(resultsTMB_list[[idx_dataset]])){
  #   plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
  #   return(NA)
  # }else{
  #   if(typeof(resultsTMB_list[[idx_dataset]]) == "character"){
  #     plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
  #     return(NA)
  #   }else{
  #     dmin1 = length(python_like_select_name(resultsTMB_list[[idx_dataset]]$par.fixed, 'beta'))/2
  #     if(full_RE){
  #       re_mat = re_vector_to_matrix(resultsTMB_list[[idx_dataset]]$par.random, dmin1)
  #       ntimes2 = length(resultsTMB_list[[idx_dataset]]$par.random)/dmin1 * 2
  #       logRmat = give_z_matrix(n_times_2 = ntimes2) %*% re_mat + 
  #         give_x_matrix(ntimes2) %*% matrix(python_like_select_name(resultsTMB_list[[idx_dataset]]$par.fixed, 'beta'), nrow=2)
  #       sim_thetas = softmax(cbind(logRmat, 0))
  #     }else{
  #       sim_thetas = softmax(cbind(sapply(1:dmin1,
  #                                         function(some_dummy_idx) give_z_matrix(length(resultsTMB_list[[idx_dataset]]$par.random) * 2) %*% resultsTMB_list[[idx_dataset]]$par.random) +
  #                                    give_x_matrix(length(resultsTMB_list[[idx_dataset]]$par.random) * 2) %*% matrix(python_like_select_name(resultsTMB_list[[idx_dataset]]$par.fixed, 'beta'), nrow=2), 0))
  #     }
  #     matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
  #     if(sum(sapply(matrices, length)) == 0){
  #       matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
  #     }
  #     ml_thetas = normalise_rw(do.call('rbind', matrices))
  # 
  #     plot(unlist(ml_thetas), unlist(sim_thetas), main=names(count_objects[idx_dataset]), cex.main=.7)
  #     abline(coef = c(0,1), col='blue', lty='dashed')
  #     return(list(ml_thetas=ml_thetas, sim_thetas=sim_thetas))
  #   }
  # }
  # function(idx_dataset, resultsTMB_list, full_RE=FALSE){
    if(is.character(idx_dataset)){
      ## if character, it's the name in count_objects
      idx_dataset_TMB = gsub("_", "", idx_dataset)
    }else{
      idx_dataset_TMB = idx_dataset
    }
    if(is.null(resultsTMB_list[[idx_dataset_TMB]])){
      plot(x=1, type = "n", main=names(count_objects[idx_dataset_TMB]), cex.main=.7)
      return(NA)
    }else{
      if(typeof(resultsTMB_list[[idx_dataset_TMB]]) == "character"){
        plot(x=1, type = "n", main=names(count_objects[idx_dataset_TMB]), cex.main=.7)
        return(NA)
      }else{
        dmin1 = length(python_like_select_name(resultsTMB_list[[idx_dataset_TMB]]$par.fixed, 'beta'))/2
        if(full_RE){
          re_mat = re_vector_to_matrix(resultsTMB_list[[idx_dataset_TMB]]$par.random, dmin1)
          ntimes2 = length(resultsTMB_list[[idx_dataset_TMB]]$par.random)/dmin1 * 2
          logRmat = give_z_matrix(n_times_2 = ntimes2) %*% re_mat + 
            give_x_matrix(ntimes2) %*% matrix(python_like_select_name(resultsTMB_list[[idx_dataset_TMB]]$par.fixed, 'beta'), nrow=2)
          sim_thetas = softmax(cbind(logRmat, 0))
        }else{
          sim_thetas = softmax(cbind(sapply(1:dmin1,
                                            function(some_dummy_idx) give_z_matrix(length(resultsTMB_list[[idx_dataset_TMB]]$par.random) * 2) %*% resultsTMB_list[[idx_dataset_TMB]]$par.random) +
                                       give_x_matrix(length(resultsTMB_list[[idx_dataset_TMB]]$par.random) * 2) %*% matrix(python_like_select_name(resultsTMB_list[[idx_dataset_TMB]]$par.fixed, 'beta'), nrow=2), 0))
        }
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        ml_thetas = normalise_rw(do.call('rbind', matrices))
        
        plot(unlist(ml_thetas), unlist(sim_thetas), main=names(count_objects[idx_dataset]), cex.main=.7)
        abline(coef = c(0,1), col='blue', lty='dashed')
        return(list(ml_thetas=ml_thetas, sim_thetas=sim_thetas))
      }
    }
}

## with only one RE
pdf("../../results/assessing_models/comparison_theta_PCAWG_M.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(1:length(results_TMB_M), function(idx_dataset){
# sapply(1:2, function(idx_dataset){
  give_plot_dataset(idx_dataset, resultsTMB_list = results_TMB_M, FALSE)
})
dev.off()

## with full RE

pdf("../../results/assessing_models/comparison_theta_PCAWG_fullRE_M.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(1:length(results_TMB_fullRE_M), function(idx_dataset){
  # sapply(1:10, function(idx_dataset){
  fullre_thetas_M = give_plot_dataset(idx_dataset, results_TMB_fullRE_M, TRUE)
})
dev.off()

## same, but without nucleotidesubstitution3
pdf("../../results/assessing_models/comparison_theta_PCAWG_fullRE_M_clean.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(names(count_objects)[which(!grepl("nucleotidesubstitution3", names(count_objects)))],
       function(idx_dataset){
  # sapply(1:2, function(idx_dataset){
  fullre_thetas_M = give_plot_dataset(idx_dataset, results_TMB_fullRE_M, TRUE)
})
dev.off()
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
# M (cont)
## Specific cases
par(mfrow=c(1,3))
special_cases_1 = sapply(c(grep('CNS-Oligonucleotidesubstitution1', names(results_TMB_M)),
                           grep('Head-SCCnucleotidesubstitution1', names(results_TMB_M)),
                           grep('Uterus-AdenoCAsignatures', names(results_TMB_M))
         ), give_plot_dataset, results_TMB_fullRE_M, TRUE)
as.numeric(melt(special_cases_1[1,1])$Var1) == melt(special_cases_1[2,1])$Var1

special_cases_2 = sapply(c(grep('Panc-Endocrinesignatures', names(results_TMB_M)),
                           grep('CNS-Oligosignatures', names(results_TMB_M)),
                           grep('CNS-PiloAstrosignatures', names(results_TMB_M))), give_plot_dataset, results_TMB_fullRE_M, TRUE)
as.numeric(melt(special_cases_2[1,1])$Var1) == melt(special_cases_2[2,1])$Var1


ggplot(cbind(ML=melt(special_cases_1[1,3]), Sim=melt(special_cases_1[2,3])),
       #data.frame(ML=unlist(special_cases_1[1,1]), Sim=unlist(special_cases_1[2,1]))
       aes(x=ML.value, y=Sim.value
           #,col=ML.Var2, alpha=0.1))+
       ))+
  geom_point()+guides(col=FALSE, alpha=FALSE)

do.call(lapply(1:3, function(i) ggplot(cbind(ML=melt(special_cases_1[1,i]),
                                     Sim=melt(special_cases_1[2,i])),
       aes(x=ML.value, y=Sim.value), labs(x='Normalised observed exposures', y='Simulated theta under model (M)'))+
  geom_point()+guides(col=FALSE, alpha=FALSE)), 'grid.arrange')

pdf("../../results/assessing_models/comparison_theta_PCAWG_M_special_cases_successful.pdf",
    height = 3, width = 10)
do.call('grid.arrange', args=list(grobs=lapply(1:3, function(i){
  ggplot(cbind(ML=melt(special_cases_1[1,i]),
               Sim=melt(special_cases_1[2,i])),
         aes(x=ML.value, y=Sim.value))+
    geom_point()+guides(col=FALSE, alpha=FALSE)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    labs(x='Normalised observed exposures', y='Simulated theta under model (M)')+
    ggtitle(names(results_TMB_M)[c(grep('CNS-Oligo', names(results_TMB_M)),
                                   grep('Eso-AdenoCA', names(results_TMB_M))[2])][i])}),
  nrow=1))
dev.off()


pdf("../../results/assessing_models/comparison_theta_PCAWG_M_special_cases_unsuccessful2.pdf",
    height = 3, width = 10)
do.call('grid.arrange', args=list(grobs=lapply(1:3, function(i){
  ggplot(cbind(ML=melt(special_cases_2[1,i]),
               Sim=melt(special_cases_2[2,i])),
         aes(x=ML.value, y=Sim.value))+
    geom_point(col='black')+guides(col=FALSE, alpha=FALSE)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", col='black')+
    labs(x='Normalised observed exposures', y='Simulated theta under model (M)')+
    ggtitle(c('Panc-Endocrinesignatures', 'CNS-Oligosignatures', 'CNS-PiloAstrosignatures')[i])
    # ggtitle(names(results_TMB_M)[c(grep('CNS-Oligo', names(results_TMB_M)),
    #                                grep('Eso-AdenoCA', names(results_TMB_M))[2])][i])
    }),
  nrow=1))
dev.off()


head(melt(special_cases_1[,1]))
head(dcast(melt(special_cases_1[,1]), Var1+Var2~L1, value.var = 'value'))
ggplot(melt(special_cases_1[,1]), aes(x=))


#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------#
# DM
pdf("../../results/assessing_models/comparison_theta_PCAWG_DM.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(1:length(results_TMB_DM), function(idx_dataset){
# sapply(6, function(idx_dataset){
  if(is.null(results_TMB_DM[[idx_dataset]])){
    plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
    return(NA)
  }else{
    if(typeof(results_TMB_DM[[idx_dataset]]) == "character"){
      plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
      return(NA)
    }else{
      if(!(results_TMB_DM[[idx_dataset]]$pdHess)){
        # no good comvergence
        plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
        return(NA)
      }else{
        beta_coefs = python_like_select_name(results_TMB_DM[[idx_dataset]]$par.fixed, 'beta')
        RE_coefs = results_TMB_DM[[idx_dataset]]$par.random
        lambda = python_like_select_name(results_TMB_DM[[idx_dataset]]$par.fixed, 'log_lambda')
        sim_thetas = replicate(20, simulate_from_DM(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
        
        ## from observed
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        ml_thetas = normalise_rw(do.call('rbind', matrices))
        ml_thetas = replicate(20, ml_thetas)
        dim(melt(sim_thetas))
        dim(melt(ml_thetas))
        plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
             main=names(count_objects[idx_dataset]), cex.main=.7)
        abline(coef = c(0,1), col='blue', lty='dashed')
      }
    }
  }
})
dev.off()

## the ones with good convergence
## !!! the overdispersion parameter has to be multiplied by the one from TMB (if using the previous implementation of lambda for DM)
results_TMB_fullRE_DM[which(unlist(sapply(results_TMB_fullRE_DM, `[`, "pdHess")) == "TRUE")]

pdf("../../results/assessing_models/comparison_theta_PCAWG_DM_fullRE.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(1:length(results_TMB_fullRE_DM), function(idx_dataset){
  # sapply(6, function(idx_dataset){
  if(is.null(results_TMB_fullRE_DM[[idx_dataset]])){
    plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
    return(NA)
  }else{
    if(typeof(results_TMB_fullRE_DM[[idx_dataset]]) == "character"){
      plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
      return(NA)
    }else{
      if(!(results_TMB_fullRE_DM[[idx_dataset]]$pdHess)){
        # no good convergence
        plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
        return(NA)
      }else{
        beta_coefs = python_like_select_name(results_TMB_fullRE_DM[[idx_dataset]]$par.fixed, 'beta')
        RE_coefs = results_TMB_fullRE_DM[[idx_dataset]]$par.random
        lambda = python_like_select_name(results_TMB_fullRE_DM[[idx_dataset]]$par.fixed, 'log_lambda')
        sim_thetas = replicate(20, simulate_from_DM_RE(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
        
        ## from observed
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        ml_thetas = normalise_rw(do.call('rbind', matrices))
        ml_thetas = replicate(20, ml_thetas)
        dim(melt(sim_thetas))
        dim(melt(ml_thetas))
        plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
             main=names(count_objects[idx_dataset]), cex.main=.7)
        abline(coef = c(0,1), col='blue', lty='dashed')
      }
    }
  }
})
dev.off()

pdf("../../results/assessing_models/comparison_theta_PCAWG_DM_diagRE.pdf")
par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
sapply(1:length(results_TMB_diagRE_DM), function(idx_dataset){
  # sapply(6, function(idx_dataset){
  if(is.null(results_TMB_diagRE_DM[[idx_dataset]])){
    plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
    return(NA)
  }else{
    if(typeof(results_TMB_diagRE_DM[[idx_dataset]]) == "character"){
      plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
      return(NA)
    }else{
      if(!(results_TMB_diagRE_DM[[idx_dataset]]$pdHess)){
        # no good convergence
        plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
        return(NA)
      }else{
        beta_coefs = python_like_select_name(results_TMB_diagRE_DM[[idx_dataset]]$par.fixed, 'beta')
        RE_coefs = results_TMB_diagRE_DM[[idx_dataset]]$par.random
        lambda = python_like_select_name(results_TMB_diagRE_DM[[idx_dataset]]$par.fixed, 'log_lambda')
        sim_thetas = replicate(20, simulate_from_DM_RE(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
        
        ## from observed
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        ml_thetas = normalise_rw(do.call('rbind', matrices))
        ml_thetas = replicate(20, ml_thetas)
        dim(melt(sim_thetas))
        dim(melt(ml_thetas))
        plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
             main=names(count_objects[idx_dataset]), cex.main=.7)
        abline(coef = c(0,1), col='blue', lty='dashed')
      }
    }
  }
})
dev.off()
#----------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------
## Simulate data with confidence intervals, to choose between M ME and DM ME

## if a CT has not run well, return NA

give_true_in_conf_int = function(idx_dataset, model, dataset_TMB){
  # sim_confint = sapply(1:20, function(idx_dataset){
  # sapply(6, function(idx_dataset){
  ## idx_dataset=9 is an example of a successful DM run
  # print(c(is.null(dataset_TMB[[idx_dataset]]),
  #         typeof(dataset_TMB[[idx_dataset]]) == "character",
  #         !(dataset_TMB[[idx_dataset]]$pdHess)))
  if(is.null(dataset_TMB[[idx_dataset]])){
    return(NA)
  }else{
    if(typeof(dataset_TMB[[idx_dataset]]) %in% c("character", "logical")){
      return(NA)
    }else{
      if(!(dataset_TMB[[idx_dataset]]$pdHess)){
        # no good convergence
        return(NA)
      }else{
        beta_coefs = python_like_select_name(dataset_TMB[[idx_dataset]]$par.fixed, 'beta')
        RE_coefs = dataset_TMB[[idx_dataset]]$par.random
        lambda = python_like_select_name(dataset_TMB[[idx_dataset]]$par.fixed, 'log_lambda')
        
        if(model == 'DM'){
          sim_thetas = replicate(1e3, simulate_from_DM_RE(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
        }else if(model == 'DM_altpar'){
          sim_thetas = replicate(1e3, simulate_from_DM_RE_altpar(beta_coefs, RE_coefs, lambda))
        }else if(model == 'M'){
          thetas_M = simulate_from_M_RE(beta_coefs, RE_coefs)
          sim_thetas = replicate(1e3, thetas_M)
        }else{
          stop('Indicate a correct <model>')
        }
        
        ## from observed
        matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
        if(sum(sapply(matrices, length)) == 0){
          matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
        }
        
        matrices = do.call('rbind', lapply(matrices, round))
        vec_observed = as.vector(matrices)
        
        mut_toll = rowSums(matrices)
        ## simulate with multinomial from these thetas
        idx_sim = 1
        multinom_draws = sapply(1:dim(sim_thetas)[3], function(simulate){
          as.vector(t(sapply(1:nrow(sim_thetas[,,idx_sim]), function(i) rmultinom( n = 1, size = mut_toll[i], prob = sim_thetas[i,,idx_sim]) )))
        })
        
        ## check if real values fall in the confidence interval, for each of the thetas
        
        multinom_draws[[1]]
        multinom_draws[[2]]
        in_conf_int = sapply(1:dim(multinom_draws)[1], function(i){
          confint_bounds = quantile(multinom_draws[i,], probs = c(0.025, 0.975))
          (vec_observed[i] >= confint_bounds[1]) & (vec_observed[i] <= confint_bounds[2])
        })
        
        ## plotting with subsample
        subset_pars = 1:length(vec_observed) #sample(1:length(vec_observed), size = 10)
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)])
        plot(ccccc$Var1, ccccc$value)
        points(1:length(subset_pars), vec_observed[subset_pars], col='red')

        ## the ones that worked
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][in_conf_int,])
        plot(ccccc$Var1, ccccc$value)
        points(1:sum(in_conf_int), vec_observed[subset_pars][in_conf_int], col='red')
        
        ## the ones that didn't work
        ccccc=melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][!in_conf_int,])
        plot(ccccc$Var1, ccccc$value)
        points(1:sum(!in_conf_int), vec_observed[subset_pars][!in_conf_int], col='red')
        
        return(in_conf_int)
        
        # ml_thetas = normalise_rw(do.call('rbind', matrices))
        # ml_thetas = replicate(20, ml_thetas)
        # dim(melt(sim_thetas))
        # dim(melt(ml_thetas))
        # plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
        #      main=names(count_objects[idx_dataset]), cex.main=.7)
        # abline(coef = c(0,1), col='blue', lty='dashed')
      }
    }
  }
}
sim_confint_M = sapply(1:length(results_TMB_fullRE_M), give_true_in_conf_int, model="M", dataset_TMB=results_TMB_fullRE_M)
names(sim_confint_M) = names(results_TMB_fullRE_M)
sim_confint_DM = sapply(1:length(results_TMB_fullRE_DM), give_true_in_conf_int, model="DM", dataset_TMB=results_TMB_fullRE_DM)
names(sim_confint_DM) = names(results_TMB_fullRE_DM)

which(grepl('Kidney-RCC.clearcell', names(results_TMB_fullRE_DM)))
names(results_TMB_fullRE_DM)[55]


all(na.omit(names(sim_confint_M) == names(sim_confint_DM))) ## check if same ct/feature combinations in all models
sim_confint_M_prop = sapply(sim_confint_M, function(i) sum(i/length(i)))
sim_confint_DM_prop = sapply(sim_confint_DM, function(i) sum(i/length(i)))
sim_confint_M_prop[is.na(sim_confint_M_prop)] = -0.2
sim_confint_DM_prop[is.na(sim_confint_DM_prop)] = -0.2

plot(sim_confint_M_prop, sim_confint_DM_prop)
abline(h = 0, lty='dashed')
abline(v = 0, lty='dashed')
abline(coef = c(0,1), lty='dashed')

confint_df = cbind.data.frame(sim_confint_M_prop,sim_confint_DM_prop, names=names(sim_confint_M_prop))
confint_df$CT = sapply(confint_df$names, function(i) gsub("nucleotidesubstitution1|signatures", "", i) )
confint_df$type = sapply(1:nrow(confint_df), function(i) gsub(confint_df$CT[i], "", confint_df$names[i]) )
ggplot(melt(confint_df, measure.vars = c('sim_confint_M_prop', 'sim_confint_DM_prop')),
       aes(x=interaction(variable, type, CT), col=variable, y=value, label=names))+geom_bar(stat='identity')+
  ggtitle("Fraction of observed values of simulations under the model with inferred parameters for the M and DM models")

ggplot(confint_df,
       aes(x=sim_confint_M_prop, y=sim_confint_DM_prop, label=CT, col=type))+
  geom_point()+geom_abline(slope = 1, intercept = 0)+geom_label_repel()+
  ggtitle("Fraction of observed values of simulations under the model with inferred parameters for the M and DM models")
ggsave("../../results/assessing_models/fraction_in_confidenceinterval_M_DM_altpar.pdf", width = 7)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------#
# LNM
simulate_from_LNM = function(beta_coefs, RE_coefs, sigma){
  mu = softmax(cbind(sapply(1:(length(beta_coefs)/2),
                                  function(some_dummy_idx){
                                    give_z_matrix(length(RE_coefs) * 2) %*% RE_coefs}) +
                             give_x_matrix(length(RE_coefs) * 2) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*exp(lambda)
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}
# simulate_from_DM = function(beta_coefs, RE_coefs, x, z){


stop('LNM is being re-computed with cov coefficients as parameters')
# pdf("../../results/assessing_models/comparison_theta_PCAWG_LNM.pdf")
# par(mfrow=c(5,5), mar=c(1.8,1.8,1.8,1.8))
# sapply(1:length(results_TMB_LNM), function(idx_dataset){
#   # sapply(6, function(idx_dataset){
#   if(is.null(results_TMB_LNM[[idx_dataset]])){
#     plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
#     return(NA)
#   }else{
#     if(typeof(results_TMB_LNM[[idx_dataset]]) == "character"){
#       plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
#       return(NA)
#     }else{
#       if(!(results_TMB_LNM[[idx_dataset]]$pdHess)){
#         # no good comvergence
#         plot(x=1, type = "n", main=names(count_objects[idx_dataset]), cex.main=.7)
#         return(NA)
#       }else{
#         beta_coefs = python_like_select_name(results_TMB_LNM[[idx_dataset]]$par.fixed, 'beta')
#         RE_coefs = results_TMB_LNM[[idx_dataset]]$par.random
#         sigma = python_like_select_name(results_TMB_LNM[[idx_dataset]]$par.fixed, 'sigma')
#         sim_thetas = replicate(20, simulate_from_DM(beta_coefs, RE_coefs, lambda))
#         
#         ## from observed
#         matrices = slot(count_objects[[idx_dataset]], 'count_matrices_active')
#         if(sum(sapply(matrices, length)) == 0){
#           matrices = slot(count_objects[[idx_dataset]], 'count_matrices_all')
#         }
#         ml_thetas = normalise_rw(do.call('rbind', matrices))
#         ml_thetas = replicate(20, ml_thetas)
#         dim(melt(sim_thetas))
#         dim(melt(ml_thetas))
#         plot(cbind(sim=melt(sim_thetas), ml=melt(ml_thetas))[,c('sim.value', 'ml.value')],
#              main=names(count_objects[idx_dataset]), cex.main=.7)
#         abline(coef = c(0,1), col='blue', lty='dashed')
#       }
#     }
#   }
# })
# dev.off()
#----------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------
## Feature correlations under the model
## positive correlations between features are possible
par(mfrow=c(1,3))
replicate(3, pairs(give_plot_dataset(1, results_TMB_fullRE_M, TRUE)$sim_thetas))
#----------------------------------------------------------------------------------------------------------

theta_observed = lapply(count_objects, function(i){
  if(is.na(i)){
    return(NA)
  }else{
    matrices = slot(i, 'count_matrices_active')
    if(sum(sapply(matrices, length)) == 0){
      matrices = slot(i, 'count_matrices_all')
    }
    ml_thetas = normalise_rw(do.call('rbind', matrices))
    ml_thetas
  }
})

correlations_observed = lapply(theta_observed, function(k){
  if(is.na(k)){
    return(NA)
  }else{
    res = outer(1:ncol(k), 1:ncol(k),
          Vectorize(function(i,j) cor(k[,i],k[,j])))
    return(res[upper.tri(res)])
  }
})

correlations_observed_patientwise = lapply(theta_observed, function(k){
  if(is.na(k)){
    return(NA)
  }else{
    res = outer(1:row(k), 1:nrow(k),
                Vectorize(function(i,j) cor(k[i,],k[j,])))
    return(res[upper.tri(res)])
  }
})

correlations_fun = function(resultsTMB_list, ...){
  lapply(1:length(resultsTMB_list), function(idx_dataset){
    if(length(resultsTMB_list[[idx_dataset]]) == 1){
      if(is.na(resultsTMB_list[[idx_dataset]])){
        return(NA)
      }else{
        return(NA)
    }
    }else{
      a = give_plot_dataset(idx_dataset, resultsTMB_list, ...)
      if(is.na(a)){
        return(NA)
      }else{
        a = a$sim_thetas
        res = outer(1:ncol(a), 1:ncol(a), Vectorize(function(i,j) cor(a[,i],a[,j])))
        return(res[upper.tri(res)])
      }
    }
  })
}

correlations_fun_patients = function(resultsTMB_list, ...){
  lapply(1:length(resultsTMB_list), function(idx_dataset){
    if(length(resultsTMB_list[[idx_dataset]]) == 1){
      if(is.na(resultsTMB_list[[idx_dataset]])){
        return(NA)
      }else{
        return(NA)
      }
    }else{
      a = give_plot_dataset(idx_dataset, resultsTMB_list, ...)
      if(is.na(a)){
        return(NA)
      }else{
        a = a$sim_thetas
        res = outer(1:nrow(a), 1:nrow(a), Vectorize(function(i,j) cor(a[i,],a[j,])))
        return(res[upper.tri(res)])
      }
    }
  })
}

## correlations between features
correlations_fullRE_M = correlations_fun(results_TMB_fullRE_M, full_RE=TRUE)
correlations_fullRE_DM = correlations_fun(results_TMB_fullRE_DM, full_RE=TRUE)
correlations_M = correlations_fun(results_TMB_M, full_RE=FALSE)
correlations_DM = correlations_fun(results_TMB_DM, full_RE=FALSE)
names(correlations_fullRE_M) = names(correlations_fullRE_DM) = names(correlations_M) = names(correlations_DM) = names(count_objects)

## do the same for patients
correlations_fullRE_M_patientwise = correlations_fun_patients(results_TMB_fullRE_M, full_RE=TRUE)
correlations_fullRE_DM_patientwise = correlations_fun(results_TMB_fullRE_DM, full_RE=TRUE)
correlations_M_patientwise = correlations_fun(results_TMB_M, full_RE=FALSE)
correlations_DM_patientwise = correlations_fun(results_TMB_DM, full_RE=FALSE)

cor1 = ggplot(melt(list(fullRE_M=correlations_fullRE_M, observed=correlations_observed)),
              aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in full RE M and observed")+theme(legend.position = "bottom")
cor2 = ggplot(melt(list(fullRE_DM=correlations_fullRE_DM, observed=correlations_observed)),
              aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in full RE DM and observed")+theme(legend.position = "bottom")
cor3 = ggplot(melt(list(M=correlations_M, observed=correlations_observed)),
              aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in M and observed")+theme(legend.position = "bottom")
cor4 = ggplot(melt(list(DM=correlations_DM, observed=correlations_observed)),
              aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in DM and observed")+theme(legend.position = "bottom")

pdf("../../results/assessing_models/correlation_features_density_model_comparison.pdf")
grid.arrange(cor1, cor2, cor3, cor4)
dev.off()

plts_features = lapply(c('signatures', 'nucleotidesubstitution1', 'nucleotidesubstitution1'), function(feature){
  plt1 = ggplot(melt(list(fullRE_M=correlations_fullRE_M[grep(feature, names(correlations_observed))],
                          observed=correlations_observed[grep(feature, names(correlations_observed))])),
         aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in fullRE M and observed")+theme(legend.position = "bottom")
  plt2 = ggplot(melt(list(fullRE_DM=correlations_fullRE_DM[grep(feature, names(correlations_observed))],
                          observed=correlations_observed[grep(feature, names(correlations_observed))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in fullRE DM and observed")+theme(legend.position = "bottom")
  plt3 = ggplot(melt(list(M=correlations_M[grep(feature, names(correlations_observed))],
                          observed=correlations_observed[grep(feature, names(correlations_observed))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in M and observed")+theme(legend.position = "bottom")
  plt4 = ggplot(melt(list(DM=correlations_DM[grep(feature, names(correlations_observed))],
                          observed=correlations_observed[grep(feature, names(correlations_observed))])),
         aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of feature correlations\n in DM and observed")+theme(legend.position = "bottom")
  list(plt1, plt2, plt3, plt4)
})

pdf("../../results/assessing_models/correlation_features_density_model_comparison_signatures.pdf",
    width = 15, height = 5)
do.call(grid.arrange, list(grobs=plts_features[[1]], nrow=1))
dev.off()

pdf("../../results/assessing_models/correlation_features_density_model_comparison_nucleotide1.pdf",
    width = 15, height = 5)
do.call(grid.arrange, list(grobs=plts_features[[2]], nrow=1))
dev.off()

pdf("../../results/assessing_models/correlation_features_density_model_comparison_nucleotide3.pdf",
    width = 15, height = 5)
do.call(grid.arrange, list(grobs=plts_features[[3]], nrow=1))
dev.off()

ggplot(melt(correlations_observed), aes(x=value, col=L1))+geom_density()+guides(col=FALSE)
head(list(correlations_fullRE_M, correlations_observed))


## Comparing the correlations between patients
plts_features_patientwise = lapply(c('signatures', 'nucleotidesubstitution1'), function(feature){
  plt1 = ggplot(melt(list(fullRE_M=correlations_fullRE_M_patientwise[grep(feature, names(correlations_observed_patientwise))],
                          observed=correlations_observed_patientwise[grep(feature, names(correlations_observed_patientwise))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of patient correlations\n in fullRE M and observed")+theme(legend.position = "bottom")
  plt2 = ggplot(melt(list(fullRE_DM=correlations_fullRE_DM_patientwise[grep(feature, names(correlations_observed_patientwise))],
                          observed=correlations_observed_patientwise[grep(feature, names(correlations_observed_patientwise))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of patient correlations\n in fullRE DM and observed")+theme(legend.position = "bottom")
  plt3 = ggplot(melt(list(M=correlations_M_patientwise[grep(feature, names(correlations_observed_patientwise))],
                          observed=correlations_observed_patientwise[grep(feature, names(correlations_observed_patientwise))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of patient correlations\n in M and observed")+theme(legend.position = "bottom")
  plt4 = ggplot(melt(list(DM=correlations_DM_patientwise[grep(feature, names(correlations_observed_patientwise))],
                          observed=correlations_observed_patientwise[grep(feature, names(correlations_observed_patientwise))])),
                aes(x=value, col=factor(L1), group=factor(L1)))+geom_density()+ggtitle("Density of patient correlations\n in DM and observed")+theme(legend.position = "bottom")
  list(plt1, plt2, plt3, plt4)
})

pdf("../../results/assessing_models/correlation_patients_density_model_comparison_signatures.pdf",
    width = 15, height = 5)
do.call(grid.arrange, list(grobs=plts_features_patientwise[[1]], nrow=1))
dev.off()


pdf("../../results/assessing_models/correlation_patients_density_model_comparison_nucleotide1.pdf",
    width = 15, height = 5)
do.call(grid.arrange, list(grobs=plts_features_patientwise[[2]], nrow=1))
dev.off()

