#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(ggplot2)
require(R.utils)
require(dplyr)
library(parallel)
library(RColorBrewer)
library(reshape2)
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
# set.seed(1234)
re_run_inference = FALSE ## use cache or not
give_summary_runs = FALSE ## whether to run the section to see what has converged, what hasn't, etc.
folder_robjs = "../../data/robjects_cache/tmb_results/"
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
TMB::compile("mm_multinomial/ME_LNM.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_LNM"))
TMB::compile("mm_multinomial/ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_multinomial"))
TMB::compile("mm_multinomial/ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_dirichletmultinomial"))
#-------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------#
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                        strsplit, split = "_")))
colnames(samples_files) = c('CT', 'type')
table(samples_files[,1], samples_files[,2])
ct = "Breast-DCIS" #samples_files[1,1]
typedata = "signatures" #samples_files[1,2]

samples_files2 = samples_files %>% filter(type != "nucleotidesubstitution3")
rownames(samples_files2) = rownames(samples_files)[samples_files$type != "nucleotidesubstitution3"]
#-------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------#
## run at random
if(re_run_inference){

  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("LNM_", "", list.files("../../data/robjects_cache/tmb_results/"))))))),
           function(idx){
  # mclapply(1:nrow(samples_files), function(idx){
    i = samples_files2[idx,]
    x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "LNM"),
                    timeout = 300, onTimeout = "warning")
    saveRDS(object = x, file=paste0("../../data/robjects_cache/tmb_results/", "LNM_", rownames(i), ".RDS"))
  })
  
  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("M_", "", list.files("../../data/robjects_cache/tmb_results/"))))))),
   function(idx){
     i = samples_files2[idx,]
     x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "M"),
                     timeout = 300, onTimeout = "warning")
     saveRDS(object = x, file=paste0("../../data/robjects_cache/tmb_results/", "M_", rownames(i), ".RDS"))
  })
  
  mclapply(sample(which(is.na(match(rownames(samples_files2),
                                    gsub(".RDS", "", gsub("DM_", "", list.files("../../data/robjects_cache/tmb_results/"))))))),
           function(idx){
             i = samples_files2[idx,]
             x = withTimeout(wrapper_run_TMB(i[1,1], i[1,2], model = "DM"),
                             timeout = 300, onTimeout = "warning")
             saveRDS(object = x, file=paste0("../../data/robjects_cache/tmb_results/", "DM_", rownames(i), ".RDS"))
           })
}

#----------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------#
# re-run those that had NaN due to bad initial values
ct = "Uterus-AdenoCA"
type= "signatures" #"nucleotidesubstitution1" #"signatures" # #
model = "DM"
rm(x)
x = wrapper_run_TMB(ct, type, model)
x
saveRDS(object = x, file=paste0("../../data/robjects_cache/tmb_results/", model, "_", ct, "_", type, ".RDS"))
#------------------------------------------------------------------------------------------------------------------------#

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
#----------------------------------------------------------------------------------------------------#

if(give_summary_runs){
  #----------------------------------------------------------------------------------------------------#
  ## checking how many errors were there
  ## a) timeout
  ## b) non-positive-definite hessian
  ## c) good convergence
  
  give_summary_of_runs(results_TMB_LNM, long_return = FALSE)
  give_summary_of_runs(results_TMB_M, long_return = FALSE)
  give_summary_of_runs(results_TMB_DM, long_return = FALSE)
  
  # long and short version of summaries
  sapply(give_summary_of_runs(results_TMB_LNM, long_return = TRUE), length)
  give_summary_of_runs(results_TMB_M, long_return = TRUE)
  
  results_TMB_M[['Biliary-AdenoCAnucleotidesubstitution1']]
  results_TMB_DM[['Biliary-AdenoCAnucleotidesubstitution1']]
  results_TMB_LNM[['Biliary-AdenoCAnucleotidesubstitution1']]
  
  ## for creating spreadsheet
  for(i in sort(unique(samples_files$CT))){cat(i, '\n')}
  
  ## see if nuc1 is positive semi-definite or has some problem
  for(i in sort(unique(samples_files$CT))){cat(i, '\n')}
  
  length(unique(samples_files$CT))
  length(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))])
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
         give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_M[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_DM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  #------------------------------------------------------------------------------------------------------------------------#
  
  #------------------------------------------------------------------------------------------------------------------------#
  for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution1", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "nucleotidesubstitution3", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  
  for(i in sapply(results_TMB_LNM[sapply(as.character(unique(samples_files$CT)), function(i) paste0(i, "", "signatures", collapse=""))],
                  give_summary_per_sample)){cat(i,'\n')}
  #------------------------------------------------------------------------------------------------------------------------#

  
  #------------------------------------------------------------------------------------------------------------------------#
  stan_results = read.table("../inference_diagnostics20200828.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
  
  for(i in get_summary_stan("M", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("M", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("M", "signatures")){cat(i,'\n')}
  
  for(i in get_summary_stan("DM", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("DM", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("DM", "signatures")){cat(i,'\n')}
  
  for(i in get_summary_stan("LNM", "nucleotidesubstitution1")){cat(i,'\n')}
  for(i in get_summary_stan("LNM", "nucleotidesubstitution3")){cat(i,'\n')}
  for(i in get_summary_stan("LNM", "signatures")){cat(i,'\n')}

  ## showing that I was not doing the parallelising correctly
  plot(stan_results$total.time..s., stan_results$real.time..s.)
  plot(density(na.omit(stan_results$real.time..s./stan_results$total.time..s.)))

  ## see if there is one run for each combination, and no more
  sort(table(apply(cbind(stan_results$CT, stan_results$features, stan_results$model), 1, paste0, collapse="")))
  ## there should be only one of each
}

#------------------------------------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------------------------------------#

betasM = load_posteriors("../../data/robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e.Rdata")
betasDM = load_posteriors("../../data/robjects_cache/betas316eb9a5-31f9-4d4b-be88-1b0e5c184286_DM.Rdata")

# some example
betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`
results_TMB_DM$`Biliary-AdenoCAsignatures`
set.seed(182)
#------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------#
## they should have the same dimensions
length(results_TMB_DM$`Biliary-AdenoCAsignatures`$par.fixed[grep("beta", names(results_TMB_DM$`Biliary-AdenoCAsignatures`$par.fixed))])
dim(betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`[,python_like_select(colnames(betasDM$posteriors_betas$`Biliary-AdenoCA_signatures_20000_DMROO`), "beta\\[2,")])

## plotting the histogram of the posteriors and the point MLE
replot = FALSE
if(replot){
  for(nme in names(results_TMB_DM)[grepl('signatures', names(results_TMB_DM))]){
    if(typeof(results_TMB_DM[[nme]]) == "character"){next}
    results_TMB_DM[[nme]]$par.fixed[grep("beta", names(results_TMB_DM[[nme]]$par.fixed))]
    
    idx_posteriors = which(grepl(nme, sapply(names(betasDM$posteriors_betas), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = ""))))
    if(! (length(idx_posteriors) == 1)){next}
    
    d_min_1 = sum(grepl('beta\\[2,', colnames(betasDM$posteriors_betas[[idx_posteriors]])))
    
    pdf(paste0("../../results/TMB/", nme, ".pdf"))
    par(mfrow=c(4,4))
    for(idx in 1:d_min_1){
      beta_val = select_slope_2(results_TMB_DM[[nme]]$par.fixed[grep("beta", names(results_TMB_DM[[nme]]$par.fixed))])[idx]
      plot(density(betasDM$posteriors_betas[[idx_posteriors]][,python_like_select(colnames(betasDM$posteriors_betas[[idx_posteriors]]), "beta\\[2,")][,idx]),
           main=paste0(idx, ' ', round(beta_val, digits = 2)))
      abline(v=beta_val,
             col='blue', lty="dashed")
    }
    dev.off()
  }
}
#------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------#
## match stan and TMB
subset_results_TMB_M = results_TMB_M[match(sapply(names(betasM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
                                           names(results_TMB_M))]
subset_results_TMB_DM = results_TMB_DM[match(sapply(names(betasDM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
      names(results_TMB_DM))]

subset_results_TMB_DM_dep = results_TMB_DM_dep[match(sapply(names(betasDM$posteriors_betas_slope), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = "")),
                                             names(results_TMB_DM_dep))]


## compare betas (note that some gave NaN for their standard error)
unlisted_TMB_list = function(TMB_list){
  sapply(1:length(TMB_list), function(dataset){
  if(typeof(TMB_list[[dataset]]) != "character"){
    TMB_list[[dataset]]$par.fixed[grepl('beta', names(TMB_list[[dataset]]$par.fixed))]
  }else{
    NA
  }
})
}

betas_TMB_M_subset = unlisted_TMB_list(subset_results_TMB_M)
unlisted_TMB_M = sapply(betas_TMB_M, select_slope_2, verbatim=FALSE)
betas_TMB_DM_subset = unlisted_TMB_list(subset_results_TMB_DM)
unlisted_TMB_DM = sapply(betas_TMB_DM_subset, select_slope_2, verbatim=FALSE)
names(unlisted_TMB_M) = subset_results_TMB_M
names(unlisted_TMB_DM) = subset_results_TMB_DM

unlisted_stan_mean_M = sapply(betasM$posteriors_betas_slope, function(i) colMeans(i))
unlisted_stan_mean_DM = sapply(betasDM$posteriors_betas_slope, function(i) colMeans(i))

unlisted_RE_coef_stan_mean_M = 0 ## I would need to get those from the cluster


## make sure they are the same length
length(unlisted_TMB_DM)
length(unlisted_stan_mean_DM)

## make sure they are the same length
length(unlisted_TMB_DM %>% unlist)
length(unlisted_stan_mean_DM %>% unlist)

## TRUE, FALSE, colour scheme
myColors <- c('#5ad3b7', '#d35a6c')
names(myColors) <- c(T, F)
colScale <- scale_colour_manual(name = "grp",values = myColors)

par(mfrow=c(1,1))
## scatterplot for multinomial draws of the mean beta from the posterior and the MLE for beta
plot(unlisted_TMB_M %>% unlist, unlisted_stan_mean_M %>% unlist,
     col=c('blue', 'red')[factor(sapply(subset_results_TMB_M, give_summary_per_sample) == "Good", levels=c(T,F))], pch=8, cex=0.5)
df_betas_slope_M = cbind.data.frame( TMB_estimate=(unlisted_TMB_M %>% unlist),
                                     stan_mean=(unlisted_stan_mean_M %>% unlist),
                                     Bool_good_convergence=factor(rep(sapply(subset_results_TMB_M, give_summary_per_sample) == "Good",
                                                                      as.vector(sapply(subset_results_TMB_M,
                                                                                       function(i) sum(names(i$par.fixed) == "beta")/2)))), ## /2 because of intercept+slope
                                     ct = rep(sapply(names(unlisted_stan_mean_M), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = " ")),
                                              sapply(unlisted_TMB_M, length)))
ggplot(df_betas_slope_M,
       aes(x=TMB_estimate, y=stan_mean, col=Bool_good_convergence))+
  geom_point()+ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  geom_abline(slope = 1, intercept = 0)+theme(legend.position = "bottom")+colScale
ggsave("../../results/TMB/betas_stan_TMB_M.pdf")

## same plot, but trying to see where it went well and where it did not
ggplot(df_betas_slope_M,
       aes(x=TMB_estimate, y=stan_mean, col=Bool_good_convergence))+
  geom_point()+ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  geom_abline(slope = 1, intercept = 0)+theme(legend.position = "bottom")+facet_wrap(.~factor(ct))+colScale
ggsave("../../results/TMB/betas_stan_TMB_M_per_ct.pdf", width = 12)


## scatterplot for dirichlet-multinomial draws of the mean beta from the posterior and the MLE for beta
plot(unlisted_TMB_DM %>% unlist, unlisted_stan_mean_DM %>% unlist, col=c('blue', 'red')[factor(sapply(subset_results_TMB_DM, give_summary_per_sample) == "Good", levels=c(T,F))], pch=8)
df_betas_slope_DM = cbind.data.frame( TMB_estimate=(unlisted_TMB_DM %>% unlist),
                                      stan_mean=(unlisted_stan_mean_DM %>% unlist),
                                      TMB_estimate_dep=(unlisted_TMB_DM_dep %>% unlist),
                                      Bool_good_convergence=factor(rep(sapply(subset_results_TMB_DM, give_summary_per_sample) == "Good",
                                                                       as.vector(sapply(subset_results_TMB_DM, function(i) sum(names(i$par.fixed) == "beta")/2)))), ## /2 because of the intercept and the slope
                                      ct = rep(sapply(names(unlisted_stan_mean_DM), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = " ")),
                                               sapply(unlisted_TMB_DM, length)))
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean, col=factor(Bool_good_convergence, levels=c(F, T))))+
  geom_point()+ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+theme(legend.position = "bottom")+
  geom_abline(intercept = 0, slope = 1)+colScale
ggsave("../../results/TMB/betas_stan_TMB_DM.pdf")

ggplot(df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE),
       aes(x=TMB_estimate, y=stan_mean, col=factor(Bool_good_convergence, levels=c(F, T))))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  theme(legend.position = "bottom")+colScale
ggsave("../../results/TMB/betas_stan_TMB_DM2.pdf")

## split by samples that had good convergence and samples that didn't
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~factor(Bool_good_convergence, levels=c(F, T)))

## checking the scatterplot of the previous version which didn't have an intercept in beta
ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of\nposterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")


## same plot, but only with good convergence
ggplot(df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE),
       aes(x=TMB_estimate, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")
ggsave("../../results/TMB/betas_stan_TMB_DM_per_ct.pdf", width = 8)

ggplot(df_betas_slope_DM,
       aes(x=TMB_estimate_dep, y=stan_mean))+
  geom_point()+geom_abline(slope = 1, intercept = 0) + ggtitle('Scatterplot of mean of posterior and TMB ML estimate')+
  theme(legend.position = "bottom")+facet_wrap(.~ct, scales = "free")

## For which cancer types is there a good correlation in DM?
df_betas_slope_DM %>% filter(Bool_good_convergence == TRUE) %>% select(ct) %>% unique() %>% unlist %>% as.vector

#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## when are random effects negligible?
RE_TMB_M = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
RE_TMB_DM = sapply(results_TMB_DM, function(i) if(typeof(i) == "list"){i$par.random}else{NA})

RE_TMB_M_subset = sapply(subset_results_TMB_M, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
RE_TMB_DM_subset = sapply(subset_results_TMB_DM, function(i) if(typeof(i) == "list"){i$par.random}else{NA})
names(RE_TMB_M_subset) = names(subset_results_TMB_M); names(RE_TMB_DM_subset) = names(subset_results_TMB_DM)

## distribution of the random effects for each cancer type, and for all DM runs
ggplot(reshape2::melt(RE_TMB_DM), aes(x=value))+geom_density()+facet_wrap(.~L1, scales = "free")

## distribution of the random effects, pooling all cancer types and patients
ggplot(reshape2::melt(RE_TMB_DM), aes(x=log(abs(value))))+geom_density()
ggplot(reshape2::melt(RE_TMB_DM), aes(x=log(abs(value))))+geom_histogram()
## large enough RE are positive (difficult to see as most RE are quite low)

## mean of the absolute random effect coefficients, for each cancer type
RE_TMB_M_mean = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){mean(abs(i$par.random))}else{NA})
RE_TMB_DM_mean = sapply(results_TMB_DM, function(i) if(typeof(i) == "list"){mean(abs(i$par.random))}else{NA})
plot(density(na.omit(RE_TMB_DM_mean)))

## how to classify cancer types between negligible and non-negligible RE?
sort(RE_TMB_DM_mean)

## is there a correlation between RE in signatures and nucleotide substitutions?
plot(vector_to_ct_list(RE_TMB_M_mean))
plot(vector_to_ct_list(RE_TMB_M_mean), ylim = c(0,1))
plot(vector_to_ct_list(RE_TMB_DM_mean))
plot(vector_to_ct_list(RE_TMB_DM_mean), ylim = c(0,0.005))
## looks like an appalling correlation

## Maybe I should instead look at the ratio of random effects to fixed effects?
RE_TMB_M_median_ratio = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){
  median(abs(i$par.random))/median(abs(select_slope_2(python_like_select_name(i$par.fixed, 'beta'), verbatim = FALSE)))}else{NA})
plot(vector_to_ct_list(RE_TMB_M_median_ratio), ylim = c(0,10)) ## this is a marginally better correlation

## compare directly the random effects variance
RE_TMB_M_REvar = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')}else{NA})
## normalised by size of fixed effects (here I am normalising by ALL beta; slope and intercept)
RE_TMB_M_REvar_norm = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')/median(abs(python_like_select_name(i$par.fixed, 'beta')))}else{NA})
RE_TMB_M_REvar_norm_slope = sapply(results_TMB_M, function(i) if(typeof(i) == "list"){python_like_select_name(i$par.fixed, 'logSigma_RE')/median(abs(select_slope_2(python_like_select_name(i$par.fixed, 'beta'), verbatim = FALSE)))}else{NA})
plot(vector_to_ct_list(RE_TMB_M_REvar)) ## this is the log
plot(vector_to_ct_list(RE_TMB_M_REvar), xlim = c(-5,4), ylim = c(-4, 3)) ## this is the log
## there doesn't seem to be a good agreement either
plot(vector_to_ct_list(RE_TMB_M_REvar_norm))
plot(vector_to_ct_list(RE_TMB_M_REvar_norm_slope))

#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## compute r*2 for betas, random effects

## we need to subset the random effects too, to match it with stan
RE_TMB_M_subset
betas_TMB_DM_subset
unlisted_stan_mean_M


#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Determine differential abundance
results_TMB_M

pvals_M = sapply(results_TMB_M, wald_TMB_wrapper)
pvals_DM = sapply(results_TMB_DM, wald_TMB_wrapper)
pvals_LNM = sapply(results_TMB_LNM, wald_TMB_wrapper)

sapply(list(pvals_M, pvals_DM, pvals_LNM), max, na.rm = TRUE)

#----------------------------------------------------------------------------------------------------------

table(pvals_M <= 0.05, pvals_DM <= 0.05, pvals_LNM <= 0.05)

## why??!
names(pvals_M)[!(names(pvals_M) %in% names(pvals_LNM))]


#----------------------------------------------------------------------------------------------------------

list_results = list(results_TMB_M, results_TMB_DM, results_TMB_LNM)
cts = sapply(list_results, function(i){
  gsub("nucleotidesubstitution1|signatures", "", names(i))
})

df_idx = lapply(1:3, function(idx){
  t(sapply(unique(cts[[idx]]),
           function(j) c(which(cts[[idx]] == j & grepl("nucleotidesubstitution1", names(list_results[[idx]]))),
                         which(cts[[idx]] == j & grepl("signature", names(list_results[[idx]]))))))
})

clnm <- colnames(df_idx[[3]])
df_idx[[3]] = do.call('rbind', df_idx[[3]])
rownames(df_idx[[3]]) = clnm


