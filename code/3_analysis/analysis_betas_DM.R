### Analyse the betas

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("../")

library(rstan)
library(reshape2)
library(gridExtra)
library(uuid)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(ggrepel)

data_inference = list.files("../data/inference/", full.names = TRUE)
data_inference = data_inference[grepl("signatures_20000_DMROO.RData", data_inference)]

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
  png(paste0("../results/betas/", nme, "_betas_DM.png"))
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


save.image(paste0("../data/robjects_cache/betas", uuid::UUIDgenerate(), "_DM.Rdata"))

ggplot(melt(cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                           nonzero_features= unlist(num_not_containing_zero),
                           features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))),
       aes(x=ct, y=value, fill=variable))+
  geom_bar(stat='identity', position = "identity", alpha=.3)
ggsave(filename = "../results/betas/all_betas_zeros_DM.png")

load=F
if(load){
  df_betas = lapply(c("../data/robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e.Rdata",
                      "../data/robjects_cache/betas316eb9a5-31f9-4d4b-be88-1b0e5c184286_DM.Rdata"),
                    function(file_to_load){
    load(file_to_load)
    df_betas_zeros = cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                                      nonzero_features= unlist(num_not_containing_zero),
                                      features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))
    df_betas_zeros_melt = melt(df_betas_zeros)
    df_betas_zeros_ratio = data.frame(ct=df_betas_zeros$ct, ratio=df_betas_zeros$nonzero_features/df_betas_zeros$features,
                                      stringsAsFactors = FALSE)
    list(df_betas_zeros=df_betas_zeros, df_betas_zeros_melt=df_betas_zeros_melt, df_betas_zeros_ratio=df_betas_zeros_ratio)
  })
  names(df_betas) = c('M', 'DM')

  ggplot(df_betas[['M']]$df_betas_zeros_melt,
         aes(x=factor(ct, levels=as.character(df_betas[['M']]$df_betas_zeros$ct[order(df_betas[['M']]$df_betas_zeros$nonzero_features)])),
             y=value, fill=variable))+labs(x="")+
    geom_bar(stat='identity', position = "identity", alpha=.3)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
  ggsave(filename = "../results/betas/all_betas_zeros.png", width = 6, height = 6)
  
  ggplot(df_betas[['DM']]$df_betas_zeros_melt,
         aes(x=factor(ct, levels=as.character(df_betas[['M']]$df_betas_zeros$ct[order(df_betas[['M']]$df_betas_zeros$nonzero_features)])),
                      y=value, fill=variable))+labs(x="")+
    geom_bar(stat='identity', position = "identity", alpha=.3)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
  ggsave(filename = "../results/betas/all_betas_zeros_DM.png", width = 6, height = 6)
  
  ggplot(df_betas[['DM']]$df_betas_zeros_ratio,
         aes(x=factor(ct, levels=as.character(ct)[order(ratio)]),
             y=ratio, label=ct ))+geom_point()+geom_label_repel()+lims(y=c(min(df_betas_zeros_ratio$ratio), 1.2))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave(filename = "../results/betas/all_betas_zeros2_DM.png", width = 14)
  
  match_zeros_ratio = cbind(DM=df_betas[['DM']]$df_betas_zeros_ratio,
                            df_betas[['M']]$df_betas_zeros_ratio[match(df_betas[['M']]$df_betas_zeros_ratio$ct,
                                                                       df_betas[['DM']]$df_betas_zeros_ratio$ct),])
  stopifnot(all(match_zeros_ratio$DM.ct == match_zeros_ratio$ct))
  ggplot(match_zeros_ratio,
         aes(x=DM.ratio,
             y=ratio, label=ct ))+geom_point(size=7, alpha=0.3, col='red')+geom_label_repel()+
    lims(y=c(min(match_zeros_ratio$ratio), 1.2))
  ggsave(filename = "../results/betas/all_betas_zeros2_both.png", width = 14)
}

as.character(df_betas_zeros_ratio$ct)[order(df_betas_zeros_ratio$ratio)]








