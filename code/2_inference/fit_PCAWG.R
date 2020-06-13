## Infence using stan_logistic_normal_multinomial_ME
## with both random effects and fixed effects

# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Sys.setenv(LANG='en')

library(rstan)
library(clusterGeneration) ## for sampling covariance matrix sigma
library(uuid)
library(optparse)
source("2_inference/helper/helper_DA_stan.R")

option_list = list(
  make_option(c("--infile"), type="character", default="", 
              help="Input ROO file", metavar="character"),
  make_option(c("--cancertype"), type="character", default="", 
              help="Cancer type, string", metavar="character"),
  make_option(c("--typedata"), type="character", default="features1", 
              help="signatures/features1/features96 [default= %default]", metavar="character"),
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

model_file_name = paste0("2_inference/stan_", opt$model, "_ME.stan")
stanc(model_file_name)

if(opt$typedata == "nucleotidesubstitution1"){
  objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
}else if(opt$typedata == "nucleotidesubstitution3"){
  objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
}else if(opt$typedata == "signatures"){
  if(is.null(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]])){
    ## no active signatures
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
  }else{
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_active")
  }
  objects_sigs_per_CT_features = lapply(objects_sigs_per_CT_features, function(i){
    rwn = rownames(i)
    .x = apply(i, 2, as.numeric)
    rownames(.x) = rwn
    round(.x)
  })
}

d = ncol(objects_sigs_per_CT_features[[1]]) ## number of signatures or features
n = nrow(objects_sigs_per_CT_features[[1]]) ## number of samples

print(n)

# covariate matrix
X = matrix(NA, nrow=2, ncol=2*n)
X[1,] = 1
X[2,] = rep(c(0,1), each=n)

Z0 = matrix(0, nrow=n, ncol=n)
diag(Z0) = 1
Z = t(rbind(Z0, Z0))

## The counts
W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])

stan_data = list(n=n,
                     d = d,
                     p=2,
                     w = W,
                     x = X,
                     Z = Z)

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 2, thin = 1)

save.image(opt$output)

