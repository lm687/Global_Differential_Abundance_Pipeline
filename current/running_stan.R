rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/")
Sys.setenv(LANG='en')

library(rstan)
library(uuid)
library(optparse)
source("2_inference/helper/helper_DA_stan.R")
source("2_inference_TMB/helper_TMB.R")

ct <- "CNS-Medullo"
opt=list()
opt$infile = paste0("../data/roo/", ct, "_signatures_ROO.RDS")
opt$output = "~/Desktop/stan_output"
opt$niterations = 2000
opt$model = "M_fullME"

fle <- opt$infile
objects_sigs_per_CT_features <- load_PCAWG(ct = ct, typedata = "signatures", simulation = F, path_to_data = "../data/")

d = ncol(objects_sigs_per_CT_features$Y)
n = nrow(objects_sigs_per_CT_features$Y)/2
X = objects_sigs_per_CT_features$x
Z = objects_sigs_per_CT_features$z
W = objects_sigs_per_CT_features$Y
model_file_name = paste0("2_inference/stan_M_fullME.stan")

rstan::stanc(file = model_file_name)

stan_data = list(n=n,
                 d = d,
                 w = W,
                 p=2,
                 x = t(X),
                 Z = t(Z))

params = c('beta', 'ularge', 'sigma_RE')
# params = c('beta', 'u', 'var_u')

# model_file_name = paste0("2_inference/stan_M_ME.stan")
# model_file_name = paste0("2_inference/stan_M_fullME.stan")

set.seed(1234)
# fit_stan <- stan(file = model_file_name, data = stan_data,
#                  iter = opt$niterations, chains = 1, cores = 4, thin = 1, pars = params,
#                  control = list(stepsize=3, adapt_delta=0.95, max_treedepth=9))
# 
# fit_stan_subset = extract(fit_stan, pars = names(fit_stan)[!grepl('ularge', names(fit_stan))])
# 
# png(paste0("~/Desktop/pairs_", ct, "_signatures_MfullRE_stan.png"),
#     height = 25, width = 25, units = "in", res = 300)
# pairs(fit_stan, pars = names(fit_stan)[!grepl('ularge', names(fit_stan))])
# dev.off()

# fit_stan_mat = as.matrix(fit_stan, tail(names(fit_stan), n=8))
# pairs(as.matrix(fit_stan_mat, tail(names(fit_stan_mat), n=8)))
# 
# pairs(fit_stan_mat, pars = names(fit_stan_mat)[!grepl('ularge', names(fit_stan_mat))])

## now fit the same, but with a DMSL, and a DMDL

## DMSL
model_file_name_DMSL = paste0("2_inference/stan_DMSL_fullME.stan")
fit_stan_DMSL <- stan(file = model_file_name_DMSL, data = stan_data,
                 iter = 25000, chains = 4, cores = 4, thin = 1, pars = params,
                 control = list(stepsize=3, adapt_delta=0.98, max_treedepth=9))
saveRDS(fit_stan_DMSL, paste0("../results/results_stan/current_fullME_DMSL_", ct, ".RDS"))
# model_file_name_DMDL = paste0("2_inference/stan_DMDL_fullME.stan")
# fit_stan_DMDL <- stan(file = model_file_name_DMSL, data = stan_data,
#                       iter = opt$niterations, chains = 1, cores = 4, thin = 1, pars = params,
#                       control = list(stepsize=3, adapt_delta=0.95, max_treedepth=9))
## http://mc-stan.org/misc/warnings.html#bfmi-low 
##' This implies that the adaptation phase of the Markov Chains did not turn out
##' well and those chains likely did not explore the posterior distribution efficiently. 

## check if: betas agree between the three models
## if DMSL is identifiable (i.e. if there are extreme correlations)
## if DMDL is identifiable (i.e. if there are extreme correlations)

# fit_stan_fullREM_betas = extract(fit_stan, pars = names(fit_stan)[grepl('beta', names(fit_stan))])
# fit_stan_fullREDMSL_betas = extract(fit_stan_DMSL, pars = names(fit_stan_DMSL)[grepl('beta', names(fit_stan_DMSL))])
# fit_stan_fullREDMDL_betas = extract(fit_stan_DMDL, pars = names(fit_stan_DMDL)[grepl('beta', names(fit_stan_DMDL))])
# 
# 
# fit_stan_fullREM_betas_mean <- sapply(fit_stan_fullREM_betas, mean)
# fit_stan_fullREDMSL_betas_mean <- sapply(fit_stan_fullREDMSL_betas, mean)
# fit_stan_fullREDMDL_betas_mean <- sapply(fit_stan_fullREDMDL_betas, mean)
# 
# pairs(cbind(fit_stan_fullREM_betas_mean, fit_stan_fullREDMSL_betas_mean, fit_stan_fullREDMDL_betas_mean))
