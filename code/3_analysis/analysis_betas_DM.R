### Analyse the betas

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("../")

## script used to be called analysis_betas_DM.R

library(rstan)
library(reshape2)
library(gridExtra)
library(uuid)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(ggrepel)

option_list = list(
  make_option(c("--model_and_feature"), type="character", default=NA, 
              help="Which model to use: M, DM, LNM, and what feature, separated by an underscore", metavar="character"),
  make_option(c("--grep_posteriors"), type="character", default=NA, 
              help="Which substring to grep to get the filenames of the relevant posteriors", metavar="character")
);

## example of grep_posteriors = signatures_20000_DMROO.RData

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_inference = list.files("../data/inference/", full.names = TRUE)
data_inference = data_inference[grepl(opt$grep_posteriors, data_inference)]

posteriors_betas = lapply(data_inference,
                    function(f){
                      if(substr(f, nchar(f), nchar(f)) == "/" | basename(f) == "NA"){
                        ## no file
                        NA
                      }else{
                        print(f)
                        load(f)
                        fit_mat = as.matrix(fit_stan)
                        if(length(fit_mat) == 0){
                          NA
                        }else{
                          tryCatch(fit_mat[,grepl('beta', colnames(fit_mat))])
                        }
		      }
                    })
names(posteriors_betas) = gsub(".RData", "", basename(data_inference))
posteriors_betas = posteriors_betas[!is.na(posteriors_betas)]

sapply(names(posteriors_betas), function(nme){
  print(nme)
  print(dim(posteriors_betas[[nme]])) 
  names_slope_betas = colnames(posteriors_betas[[nme]])[c(F,T)]
  names_intersect_betas = colnames(posteriors_betas[[nme]])[c(T,F)]
  png(paste0("../results/betas/", nme, "_betas_", opt$model_and_feature, ".png"))
  grid.arrange(bayesplot::mcmc_areas(posteriors_betas[[nme]], pars = names_slope_betas)+ggtitle('Slope'),
               bayesplot::mcmc_areas(posteriors_betas[[nme]], pars = names_intersect_betas)+ggtitle('Intersect'))
  dev.off()
})

posteriors_betas_slope = lapply(posteriors_betas, function(i) i[,c(F,T)])

num_not_containing_zero = lapply(posteriors_betas_slope, function(p){
  posteriors_slopes_quant = apply(p, 2, quantile, c(0.025, 0.975))
  posteriors_slopes_quant_bool = apply(posteriors_slopes_quant, 2, function(i){
    (i[1] < 0) & (i[2] > 0)
  } )
  dim(p)[2] - sum(posteriors_slopes_quant_bool)
})


save.image(paste0("../data/robjects_cache/betas", uuid::UUIDgenerate(), "_", opt$model_and_feature, ".Rdata"))

ggplot(melt(cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                           nonzero_features= unlist(num_not_containing_zero),
                           features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))),
       aes(x=ct, y=value, fill=variable))+
  geom_bar(stat='identity', position = "identity", alpha=.3)
ggsave(filename = "../results/betas/all_betas_zeros_", opt$model_and_feature, ".png")








