## Infence using stan_logistic_normal_multinomial_ME
## with both random effects and fixed effects

# rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# Sys.setenv(LANG='en')

library(rstan)
library(clusterGeneration) ## for sampling covariance matrix sigma
library(uuid)
library(optparse)
source("helper/helper_DA_stan.R")

option_list = list(
  make_option(c("--cancer_type"), type="character", default="", 
              help="Number of samples [default= %default]", metavar="character"),
  make_option(c("--type_data"), type="character", default="features1", 
              help="signatures/features1/features96 [default= %default]", metavar="character"),
  make_option(c("--type_split_into_groups"), type="character", default="features1",
              help="hardthreshCCF1/subclonalPCAWG/replicationtiming [default= %hardthreshCCF1]", metavar="character")
  
); 
# setwd("~/DifferentialAbundance_stan_link/")
# opt=list(); opt$cancer_type = "SKCM-US"; opt$type_data = "signatures"; opt$type_split_into_groups = "subclonalPCAWG"
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

model_file_name_DM = "stan_dirichlet_multinomial_ME_new.stan"
stanc(model_file_name_DM)

Nits = 25000

# stopifnot(opt$type_data == "features1")

source("helper/header_load_ROO.R")


d = ncol(objects_sigs_per_CT_features[[1]]) ## number of signatures or features
n = nrow(objects_sigs_per_CT_features[[1]]) ## number of samples

# covariate matrix
X = matrix(NA, nrow=2, ncol=2*n)
X[1,] = 1
X[2,] = rep(c(0,1), each=n)

Z0 = matrix(0, nrow=n, ncol=n)
diag(Z0) = 1
Z = t(rbind(Z0, Z0))

## The counts
W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])

print(head(W))

stan_data_LNM = list(n=n,
                     d = d,
                     p=2,
                     w = W,
                     x = X,
                     Z = Z)

fit_LNM <- stan(file = model_file_name_DM, data = stan_data_LNM, algorithm = "NUTS",
                warmup = 1000, iter = Nits, chains = 4, cores = 4, thin = 1, control=list(adapt_delta=0.8, max_treedepth=20))

flder_out = paste0("../../../out/DA_stan/DM_PCAWG_", Nits, "_", opt$type_split_into_groups, "/")
system(paste0("mkdir -p ", flder_out))
save.image(paste0(flder_out, "DM_PCAWG_", "Nits", Nits, "_", opt$type_data, '_', opt$cancer_type, '_',
                  uuid::UUIDgenerate(), ".Rdata"))



