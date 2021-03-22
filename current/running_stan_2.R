rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/")
Sys.setenv(LANG='en')

library(rstan)
library(uuid)
library(optparse)
source("2_inference/helper/helper_DA_stan.R")
source("2_inference_TMB/helper_TMB.R")

opt=list()
opt$infile = "../data/roo/CNS-Medullo_signatures_ROO.RDS"
opt$output = "~/Desktop/stan_output"
opt$niterations = 2000
opt$model = "M_fullME"

ct = "CNS-Medullo"

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
fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = opt$niterations, chains = 1, cores = 4, thin = 1, pars = params,
                 control = list(stepsize=3, adapt_delta=0.95, max_treedepth=9))

fit_stan_subset = extract(fit_stan, pars = names(fit_stan)[!grepl('ularge', names(fit_stan))])

png(paste0("~/Desktop/pairs_", ct, "_signatures_MfullRE_stan.png"),
    height = 25, width = 25, units = "in", res = 300)
pairs(fit_stan, pars = names(fit_stan)[!grepl('ularge', names(fit_stan))])
dev.off()

fit_stan = as.matrix(fit_stan, tail(names(fit_stan), n=8))
pairs(as.matrix(fit_stan, tail(names(fit_stan), n=8)))

pairs(fit_stan, pars = names(fit_stan)[!grepl('ularge', names(fit_stan))])
