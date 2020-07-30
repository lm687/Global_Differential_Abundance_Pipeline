## Infence using stan_logistic_normal_multinomial_ME
## with both random effects and fixed effects

# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("../")
# Sys.setenv(LANG='en')

library(rstan)
library(uuid)
library(optparse)
source("2_inference/helper/helper_DA_stan.R")

option_list = list(
  make_option(c("--infile"), type="character", default="", 
              help="Input ROO file", metavar="character"),
  make_option(c("--output"), type="character", default="", 
              help="Output file", metavar="character"),
  make_option(c("--niterations"), type="numeric", default="", 
              help="Number of iterations in MCMC", metavar="character"),
  make_option(c("--model"), type="character", default="", 
              help="Model: M, DM, LNM", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Nits = opt$niterations

fle <- opt$infile
if(file.exists(fle)){
  objects_sigs_per_CT_features <- readRDS(fle)
}else{
  cat('Object for cancer type', opt$cancertype, 'was not found')
}
d = objects_sigs_per_CT_features$d ## number of signatures or features
n = objects_sigs_per_CT_features$n ## number of samples
X = objects_sigs_per_CT_features$X
Z = objects_sigs_per_CT_features$Z
objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features$objects_counts,"count_matrices_all")
W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])

model_file_name = paste0("2_inference/stan_", opt$model, "_ME.stan")
rstan_options(auto_write = TRUE)
stanc(model_file_name)

stan_data = list(n=n,
                 d = d,
                 p=2,
                 w = W,
                 x = X,
                 Z = Z)

params = c('beta', 'u', 'sigma_u')

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 4, thin = 1, pars = params, control = list(stepsize=3, adapt_delta=0.95, max_treedepth=9))

save.image(opt$output)

