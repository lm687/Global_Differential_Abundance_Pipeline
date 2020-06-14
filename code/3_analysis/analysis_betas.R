### Analyse the betas

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("../")

library(rstan)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(uuid)
library(bayesplot)

data_inference = list.files("../data/inference/", full.names = TRUE)
data_inference = data_inference[grepl("20000ROO", data_inference)]

posteriors_betas = lapply(data_inference,
                    function(f){
                      if(substr(f, nchar(f), nchar(f)) == "/" | basename(f) == "NA"){
                        ## no file
                        NA
                      }else{
                        print(f)
                        load(f)
                        fit_mat = as.matrix(fit_stan)
                        tryCatch(fit_mat[,grepl('beta', colnames(fit_mat))])
                        # tryCatch(extract(fit_stan)$beta)
                      }
                    })
names(posteriors_betas) = gsub(".RData", "", basename(data_inference))
posteriors_betas = posteriors_betas[!is.na(posteriors_betas)]

sapply(names(posteriors_betas), function(nme){
  names_slope_betas = colnames(posteriors_betas[[1]])[c(F,T)]
  names_intersect_betas = colnames(posteriors_betas[[1]])[c(T,F)]
  png(paste0("../results/betas/", nme, "_betas.png"))
  grid.arrange(bayesplot::mcmc_areas(posteriors_betas[[1]], pars = names_slope_betas)+ggtitle('Slope'),
               bayesplot::mcmc_areas(posteriors_betas[[1]], pars = names_intersect_betas)+ggtitle('Intersect'))
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


save.image(paste0("../data/robjects_cache/betas", uuid::UUIDgenerate(), ".Rdata"))

ggplot(melt(cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                           nonzero_features= unlist(num_not_containing_zero),
                           features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))),
       aes(x=ct, y=value, fill=variable))+
  geom_bar(stat='identity', position = "identity", alpha=.3)
