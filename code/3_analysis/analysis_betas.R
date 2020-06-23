### Analyse the betas




library(rstan)
library(reshape2)
library(gridExtra)
library(uuid)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(ggrepel)

debug = FALSE
if(debug){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../")
  opt = list(); opt$files_posteriors = "../data/inference/Biliary-AdenoCA_signatures_20000_MROO.RData	../data/inference/Skin-Melanoma.acral_signatures_20000_MROO.RData ../data/inference/Kidney-RCC.papillary_signatures_20000_MROO.RData ../data/inference/Skin-Melanoma.cutaneous_signatures_20000_MROO.RData"
  opt$model = 'M'
}else{
  option_list = list(
    make_option(c("--files_posteriors"), type="character", default=NA, 
                help="File with the posterior, with directory included", metavar="character"),
    make_option(c("--cancer_type"), type="character", default=NA, 
                help="File with the posterior, with directory included", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

opt$files_posteriors = strsplit(opt$files_posteriors, " ")[[1]]
data_inference = opt$files_posteriors
model = opt$model = 'M'

# data_inference = list.files("../data/inference/", full.names = TRUE)
# data_inference = data_inference[grepl("20000_MROO", data_inference)]

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
  png(paste0("../results/betas/", nme, "_betas.png"))
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


save.image(paste0("../data/robjects_cache/betas", uuid::UUIDgenerate(), ".Rdata"))

ggplot(melt(cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                           nonzero_features= unlist(num_not_containing_zero),
                           features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))),
       aes(x=ct, y=value, fill=variable))+
  geom_bar(stat='identity', position = "identity", alpha=.3)
ggsave(filename = "../results/betas/all_betas_zeros.png")

load=F
if(load){
  load("../data/robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e.Rdata")
  posteriors_betas_M = posteriors_betas
  num_not_containing_zero_M = num_not_containing_zero
  load("../data/robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e.Rdata")
  posteriors_betas_DM = posteriors_betas
  posteriors_betas = list(posteriors_betas_M, posteriors_betas_DM)
  posteriors_betas = list(posteriors_betas_M, posteriors_betas_DM)
  num_not_containing_zero_DM = num_not_containing_zero
  df_betas_zeros = cbind.data.frame(ct= sapply( names(posteriors_betas), function(i) strsplit(i, '_')[[1]][1]),
                   nonzero_features= unlist(num_not_containing_zero),
                   features=sapply(posteriors_betas_slope, function(i) dim(i)[2]))
  df_betas_zeros_melt = melt(df_betas_zeros)
  df_betas_zeros_ratio = data.frame(ct=df_betas_zeros$ct, ratio=df_betas_zeros$nonzero_features/df_betas_zeros$features,
                                    stringsAsFactors = FALSE)
  
  ggplot(df_betas_zeros_melt,
         aes(x=ct, y=value, fill=variable))+
    geom_bar(stat='identity', position = "identity", alpha=.3)
  
  ggplot(df_betas_zeros_ratio,
         aes(x=factor(ct, levels=as.character(ct)[order(ratio)]),
             y=ratio, label=ct ))+geom_point()+geom_label_repel()+lims(y=c(min(df_betas_zeros_ratio$ratio), 1.2))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave(filename = "../results/betas/all_betas_zeros2.png", width = 14)
  
  as.character(df_betas_zeros_ratio$ct)[order(df_betas_zeros_ratio$ratio)]
  
  ## load tracksig
  tracksig = read.csv("../data/restricted/tracksig/changepoints_stats_tracksig.csv", stringsAsFactors = FALSE)
  tracksig = tracksig %>% group_by(type) %>% dplyr::summarize(count = n(), bool_changepoints=sum(n_changepoints > 0))%>%
    mutate(tracksig_frac= bool_changepoints/count  )
  tracksig = cbind(df_betas_zeros_ratio, (tracksig[match(as.character(df_betas_zeros_ratio$ct), tracksig$type),]))
  
  ggplot(tracksig, aes(x=ratio, y=tracksig_frac, label=ct))+geom_point()+geom_label_repel()+
    labs(x='Ratio of non-zero coefficients', y='Fraction of TrackSig samples with some changepoint')
  ggsave(paste0("../results/betas/betas_tracksig_comparison_", model, '.png'), width = 8, height = 8)
}
