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
# opt=list(); opt$cancer_type = "BLCA-US"; opt$type_data = "signatures"; opt$type_split_into_groups = "hardthreshCCF1"

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Nits = 25000

mode_run = as.character(read.table("~/mode_run", sep = " ", stringsAsFactors = FALSE)[3])
if(mode_run == "local"){
  flder_phd = "/Users/morril01/Documents/PhD/CDA_in_Cancer/"
}else if(mode_run == 'CRUK'){
  flder_phd = "/mnt/scratcha/fmlab/morril01/git_phd/" 
}else if(mode_run == 'HPC'){
  flder_phd = "/rds/user/lm687/hpc-work/git_phd/" 
}
# stopifnot(opt$type_data == "features1")

if(opt$type_data == "features1"){
  if(opt$type_split_into_groups == "hardthreshCCF1"){
    # savefolder_features = "/Users/morril01/Documents/PhD/CDA_in_Cancer/out/ROO_PCAWG/clonal_subclonal_features1_ROO/" ## cruk laptop
    savefolder_features = paste0(flder_phd, "/out/ROO_PCAWG/clonal_subclonal_features1_ROO/")
  }else if(opt$type_split_into_groups == "subclonalPCAWG"){
    savefolder_features = paste0(flder_phd, "out/ROO_PCAWG/deconvolution_clonal_subclonal_features1_ROO/")
  }
}else if(opt$type_data == "features3"){
  if(opt$type_split_into_groups == "hardthreshCCF1"){
    # savefolder_features = "/Users/morril01/Documents/PhD/CDA_in_Cancer/out/ROO_PCAWG/clonal_subclonal_features3_ROO/" ## cruk laptop
    savefolder_features = paste0(flder_phd, "out/ROO_PCAWG/clonal_subclonal_features3_ROO/")
  }else if(opt$type_split_into_groups == "subclonalPCAWG"){
    savefolder_features = paste0(flder_phd, "out/ROO_PCAWG/deconvolution_clonal_subclonal_features3_ROO/")
  }
}else if(opt$type_data == "signatures"){
  if(opt$type_split_into_groups == "hardthreshCCF1"){
    savefolder_features = paste0(flder_phd, "out/ROO_PCAWG/clonal_subclonal_ROO/")
  }else if(opt$type_split_into_groups == "subclonalPCAWG"){
    savefolder_features = paste0(flder_phd, "out/ROO_PCAWG/deconvolution_clonal_subclonal_ROO/")
  }
}

cat('The folder for ROO files is ', savefolder_features, '\n')

objects_sigs_per_CT_features <- list()
fle <- paste0(savefolder_features, opt$cancer_type, '_ROOSigs.RDS')
if(file.exists(fle)){
  objects_sigs_per_CT_features <- readRDS(fle)
}else{
  cat('Object for cancer type', opt$cancer_type, 'was not found')
}

model_file_name_DM = "stan_multinomial_ME.stan"
stanc(model_file_name_DM)

if(opt$type_data == "features1"){
  objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
}else if(opt$type_data == "features3"){
  objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
}else if(opt$type_data == "signatures"){
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

# covariate matrix
X = matrix(NA, nrow=2, ncol=2*n)
X[1,] = 1
X[2,] = rep(c(0,1), each=n)

Z0 = matrix(0, nrow=n, ncol=n)
diag(Z0) = 1
Z = t(rbind(Z0, Z0))

## The counts
W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])

stan_data_LNM = list(n=n,
                     d = d,
                     p=2,
                     w = W,
                     x = X,
                     Z = Z)

fit_LNM <- stan(file = model_file_name_DM, data = stan_data_LNM,
                warmup = 1000, iter = Nits, chains = 4, cores = 2, thin = 1)

flder_out = paste0("../../../out/DA_stan/M_PCAWG_", Nits, '_', opt$type_split_into_groups, "/")
system(paste0("mkdir -p ", flder_out))
save.image(paste0(flder_out, "M_PCAWG_", "Nits", Nits, "_", opt$type_data, '_', opt$cancer_type, '_',
                  uuid::UUIDgenerate(), ".Rdata"))

