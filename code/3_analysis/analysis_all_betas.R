## taken from analysis_betas_DM.R and analysis_betas.R
## Analysos once files with beta (much smaller than the posteriors) have been created

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../") ## because if I do it from the cluster, I am in code/

require(reshape2)
require(ggplot2)

## load the betas from June, prior to me deleting all runs accidentally in August, and the re-runs in August

df_betas = lapply(c("../data/robjects_cache/betas91ecb3fe-4ff0-4e91-b6f0-a2eaf027f91e_M_signatures.Rdata",
                    "../data/robjects_cache/betas3a901d5c-49b1-4c6b-a52d-fd9021cd7b2a_M_signatures.Rdata",
                    "../data/robjects_cache/betas316eb9a5-31f9-4d4b-be88-1b0e5c184286_DM_signatures.Rdata",
                    "../data/robjects_cache/betas58a02ef3-12e1-485a-a8b0-7f34a689bc4e_DM_signatures.Rdata"),
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
names(df_betas) = c('M_signatures_old', 'M_signatures_new', 'DM_signatures_old', 'DM_signatures_new')

levels_ct = c(as.character(df_betas[['M_signatures_new']]$df_betas_zeros$ct[order(df_betas[['M_signatures_new']]$df_betas_zeros$nonzero_features)]),
              as.character(df_betas[['DM_signatures_new']]$df_betas_zeros$ct[!(df_betas[['DM_signatures_new']]$df_betas_zeros$ct %in% df_betas[['M_signatures_new']]$df_betas_zeros$ct)]))

ggplot(df_betas[['M_signatures_new']]$df_betas_zeros_melt,
       aes(x=factor(ct, levels=levels_ct),
           y=value, fill=variable))+labs(x="")+
  geom_bar(stat='identity', position = "identity", alpha=.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(filename = "../results/betas/all_betas_zeros_Aug.png", width = 6, height = 6)

ggplot(df_betas[['M_signatures_old']]$df_betas_zeros_melt,
       aes(x=factor(ct, levels=as.character(df_betas[['M_signatures_old']]$df_betas_zeros$ct[order(df_betas[['M_signatures_old']]$df_betas_zeros$nonzero_features)])),
           y=value, fill=variable))+labs(x="")+
  geom_bar(stat='identity', position = "identity", alpha=.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(filename = "../results/betas/all_betas_zeros.png", width = 6, height = 6)

match_M_signatures = cbind(df_betas[['M_signatures_old']]$df_betas_zeros_ratio[match(df_betas[['M_signatures_new']]$df_betas_zeros_ratio$ct,
      df_betas[['M_signatures_old']]$df_betas_zeros_ratio$ct),],
      new=df_betas[['M_signatures_new']]$df_betas_zeros_ratio)
match_DM_signatures = cbind(df_betas[['DM_signatures_old']]$df_betas_zeros_ratio[match(df_betas[['DM_signatures_new']]$df_betas_zeros_ratio$ct,
                                                                                df_betas[['DM_signatures_old']]$df_betas_zeros_ratio$ct),],
                            new=df_betas[['DM_signatures_new']]$df_betas_zeros_ratio)

plot(match_M_signatures$ratio, match_M_signatures$ratio)
plot(match_DM_signatures$ratio, match_DM_signatures$ratio)

ggplot(df_betas[['DM_signatures_old']]$df_betas_zeros_melt,
       aes(x=factor(ct, levels=as.character(df_betas[['M_signatures_old']]$df_betas_zeros$ct[order(df_betas[['M_signatures_old']]$df_betas_zeros$nonzero_features)])),
           y=value, fill=variable))+labs(x="")+
  geom_bar(stat='identity', position = "identity", alpha=.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(filename = "../results/betas/all_betas_zeros_DM_old.png", width = 6, height = 6)

ggplot(df_betas[['DM_signatures_new']]$df_betas_zeros_melt,
       aes(x=droplevels(factor(ct, levels=levels_ct)),
           y=value, fill=variable))+labs(x="")+
  geom_bar(stat='identity', position = "identity", alpha=.3)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(filename = "../results/betas/all_betas_zeros_DM_Aug.png", width = 6, height = 6)

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

as.character(df_betas_zeros_ratio$ct)[order(df_betas_zeros_ratio$ratio)]


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