rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("recover_COSMIC_signatures.R")
source("../../2_inference_TMB/helper_TMB.R")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)

v <- 'v2'
if(v == 'v1'){
  vec_sigs <- c("SBS1",   "SBS2",   "SBS3",   "SBS5",   "SBS9",   "SBS13",  "SBS15",  "SBS17a", "SBS17b", "SBS18",
               "SBS20",  "SBS21",  "SBS26",  "SBS28",  "SBS40" , "SBS41",  "SBS43",  "SBS44",  "SBS51",  "SBS58")
  add=''
}else if(v == 'v2'){
  sigs_cosmic <- read.table(paste0( "../../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                            stringsAsFactors = FALSE, sep = ',', header = TRUE)
  vec_sigs <- colnames(sigs_cosmic)[-c(1,2)]
  
  add = '_allsigs'
}
current_exposures <- list()

for(it_sig in vec_sigs){
  initial_abundance <- c(rep(0,length(vec_sigs)))
  initial_abundance[which(vec_sigs == it_sig)] <- 1
  # current_exposures[[1]] <- MCMCpack::rdirichlet(n=200, rep(1/length(vec_sigs), length(vec_sigs)))
  current_exposures[[1]] <- MCMCpack::rdirichlet(n=200, initial_abundance)
  current_exposures[[1]] <- sweep(current_exposures[[1]], 1, 260, '*')
  
  maxit <- 100
  for(it in 2:maxit){
    cat('Iteration ', it, '\n')
    .x <- give_plot_bleeding(names_sigs=vec_sigs,
                       abundances = NA, return_dataframe = T, rel_path = "../../../", exposures_input = current_exposures[[it-1]])
    .exposuressigsQP <- matrix(as.numeric(.x$exposuressigsQP), ncol=length(vec_sigs))
    # .x <- dcast(.x[,c('exposuressigsQP', 'sig', 'sample')], sig~sample, value.var = 'exposuressigsQP')
    # .x <- .x[,-1]
    # current_exposures[[it]] <- t(.x)
    current_exposures[[it]] <- .exposuressigsQP
  }

  # par(mfrow=c(4, 5))
  # for(k in 1:length(vec_sigs)){
  #   plot(sapply(current_exposures, function(i) i[,k]), type='l', main=vec_sigs[k])
  # }
  
  # current_exposures_normalised <- lapply(current_exposures, normalise_rw)
  
  current_exposures_melt <- melt(current_exposures)
  current_exposures_melt$Var2 = vec_sigs[current_exposures_melt$Var2]
  # current_exposures_normalised_melt <- melt(current_exposures_normalised)
  
  levels_sigs <- paste0('SBS', gtools::mixedsort(gsub('SBS', '', unique(current_exposures_melt$Var2))))
  # ggplot(current_exposures_melt,
  #        aes(x=L1, y=value, group=Var1))+geom_line(alpha=0.02)+
  #   facet_wrap(.~factor(Var2, levels=levels_sigs))+
  #   theme_bw()+labs(x='Iteration', y='Exposure')
  # ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, ".pdf"),
  #        height = 8, width = 8)
  
  current_exposures_melt_grouped <- current_exposures_melt %>% group_by(Var2, L1) %>%
    summarise(mean_value=mean(value), se_value=sd(value),
              perc5=quantile(value, probs = c(0.05)), perc95=quantile(value, probs = c(0.95)))
  final_mean_exposure <- current_exposures_melt_grouped[current_exposures_melt_grouped$L1 == max(current_exposures_melt_grouped$L1),]
  top95 <- final_mean_exposure$Var2[final_mean_exposure$perc95 > 20]
  current_exposures_melt_grouped$group = paste0((current_exposures_melt_grouped$Var2 == it_sig), (current_exposures_melt_grouped$Var2 %in% top95))
  current_exposures_melt_grouped$group[grepl('^TRUE*', current_exposures_melt_grouped$group)] = it_sig
  current_exposures_melt_grouped$group[current_exposures_melt_grouped$group == 'FALSETRUE'] = 'High increase'
  current_exposures_melt_grouped$group[current_exposures_melt_grouped$group == 'FALSEFALSE'] = 'Moderate increase'
  current_exposures_melt_grouped$group <- factor(current_exposures_melt_grouped$group,
                                                 levels=c(it_sig, 'High increase', 'Moderate increase'))
  
  current_exposures_melt_grouped$Var2 = factor(current_exposures_melt_grouped$Var2,
                                               levels=final_mean_exposure$Var2[order(final_mean_exposure$mean_value, decreasing = T)])
  ggplot(current_exposures_melt_grouped, aes(x=L1, y=mean_value, group=Var2, col=(Var2)))+
    # geom_ribbon(aes(ymin=mean_value-se_value, ymax=mean_value+se_value, fill=Var2), alpha=0.2)+
    geom_ribbon(aes(ymin=perc5, ymax=perc95, fill=Var2), alpha=0.2, colour=NA)+
    geom_line()+facet_wrap(.~group, scales = "free_y")+
    # geom_label_repel(data = current_exposures_melt_grouped[current_exposures_melt_grouped$L1 == 100,],
    #            aes(x=100, y=mean_value, label=Var2), size=1.8, label.size=NA, fill=NA,
    #            direction = "y", nudge_x = 20, max.overlaps = Inf)+lims(x=c(0, 140))+
    guides(fill=F,  col=guide_legend(nrow=7,byrow=TRUE, title = ''))+theme_bw()+labs(x='Iteration', y='Exposure')+
    theme(legend.position = "bottom",  legend.spacing.x = unit(.1, 'cm'), legend.spacing.y = unit(0, 'cm'))
  ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_ribbon.pdf"),
         height = 4, width = 8)
  
  tikzDevice::tikz(filename = paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_ribbon.tex"),
                   height = 4, width = 8)
  print(ggplot(current_exposures_melt_grouped, aes(x=L1, y=mean_value, group=Var2, col=(Var2)))+
    # geom_ribbon(aes(ymin=mean_value-se_value, ymax=mean_value+se_value, fill=Var2), alpha=0.2)+
    geom_ribbon(aes(ymin=perc5, ymax=perc95, fill=Var2), alpha=0.2, colour=NA)+
    geom_line()+facet_wrap(.~group, scales = "free_y")+
    # geom_label_repel(data = current_exposures_melt_grouped[current_exposures_melt_grouped$L1 == 100,],
    #            aes(x=100, y=mean_value, label=Var2), size=1.8, label.size=NA, fill=NA,
    #            direction = "y", nudge_x = 20, max.overlaps = Inf)+lims(x=c(0, 140))+
    guides(fill=F,  col=guide_legend(nrow=7,byrow=TRUE, title = ''))+theme_bw()+labs(x='Iteration', y='Exposure')+
    theme(legend.position = "bottom",  legend.spacing.x = unit(.1, 'cm'), legend.spacing.y = unit(0, 'cm')))
  dev.off()
  
  tikzDevice::tikz(filename = paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_ribbon2.tex"),
                   height = 2.5, width = 6)
  print(ggplot(current_exposures_melt_grouped, aes(x=L1, y=mean_value, group=Var2, col=(Var2)))+
    # geom_ribbon(aes(ymin=mean_value-se_value, ymax=mean_value+se_value, fill=Var2), alpha=0.2)+
    geom_ribbon(aes(ymin=perc5, ymax=perc95, fill=Var2), alpha=0.2, colour=NA)+
    geom_line()+facet_wrap(.~group, scales = "free_y")+
    # geom_label_repel(data = current_exposures_melt_grouped[current_exposures_melt_grouped$L1 == 100,],
    #            aes(x=100, y=mean_value, label=Var2), size=1.8, label.size=NA, fill=NA,
    #            direction = "y", nudge_x = 20, max.overlaps = Inf)+lims(x=c(0, 140))+
    guides(fill=F,  col=F)+theme_bw()+labs(x='Iteration', y='Exposure')+
    theme(legend.position = "bottom",  legend.spacing.x = unit(.1, 'cm'), legend.spacing.y = unit(0, 'cm')))
  dev.off()
  
  ggplot(current_exposures_melt_grouped, aes(x=L1, y=mean_value, group=Var2, col=(Var2)))+
          geom_ribbon(aes(ymin=perc5, ymax=perc95, fill=Var2), alpha=0.2, colour=NA)+
          geom_line()+facet_wrap(.~group, scales = "free_y")+
    geom_label_repel(data = current_exposures_melt_grouped[current_exposures_melt_grouped$L1 == 100,],
               aes(x=100, y=mean_value, label=ifelse(perc95>40, as.character(Var2), NA)), label.size=NA, fill=NA,
               direction = "y", nudge_x = 10, max.overlaps = Inf, size=3)+lims(x=c(0, 120))+
    guides(fill=F,  col=F)+theme_bw()+labs(x='Iteration', y='Exposure')+
          theme(legend.position = "bottom",  legend.spacing.x = unit(.1, 'cm'), legend.spacing.y = unit(0, 'cm'))
  ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_ribbon2.png"),
         height = 2.5, width = 6)
  
  fill_with_one <- function(which_one, n=length(vec_sigs)){
    .x <- rep(0, n)
    .x[which_one] <- 1
    .x
  }
  vec_sigs[sapply(1:length(current_exposures), function(i) which.max(colSums(current_exposures[[i]])))]
  which_max_df = melt(sapply(1:length(current_exposures), function(i) fill_with_one(which.max(colSums(current_exposures[[i]])))))
  which_max_df$Var1 = vec_sigs[which_max_df$Var1]
  which_max_df <- which_max_df[which_max_df$value == 1,]
  ggplot(which_max_df,
         aes(y=factor(Var1, gtools::mixedsort(unique(Var1))),
             # x=factor(Var2, levels=sort(unique(Var2), decreasing = F)),
             x=Var2,
             fill=value, group=1))+
    geom_tile()+geom_path()+theme_bw()+
    # theme(axis.text.x=element_text(angle = 45, hjust = 1))+
    labs(y='Signatures', x='Iteration')+guides(fill=F)
  ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_trajectory.pdf"),
         height = 2.5, width = 6)
  
  per_patient_trajectory <- F
  if(per_patient_trajectory){
    which_max_df_persample = melt(lapply(1:length(current_exposures),
                                         function(i) apply(current_exposures[[i]], 1, function(p) fill_with_one(which.max(p)))))
    which_max_df_persample$Var1 = vec_sigs[which_max_df_persample$Var1]
    which_max_df_persample <- which_max_df_persample[which_max_df_persample$value == 1,]
    ggplot(which_max_df_persample,
           aes(y=factor(Var1, gtools::mixedsort(unique(Var1))),
               x=L1,
               fill=value, group=1))+
      geom_tile(alpha=0.2)+
      geom_path(alpha=0.2)+theme_bw()+
      labs(y='Signatures', x='Iteration')+guides(fill=F)
    
  }
  ## NOT YET RUN FOR ALL SAMPLES
  
  AUC_mat <-  Reduce('+', current_exposures)
  AUC_signature <- sum(AUC_mat[,which(vec_sigs == it_sig)])/sum(AUC_mat)
  AUC_signatures <- sapply(vec_sigs, function(it_sig) sum(AUC_mat[,which(vec_sigs == it_sig)])/sum(AUC_mat))
  
  ggplot(data.frame(AUC_signatures=AUC_signatures,
                    sig=names(AUC_signatures))[AUC_signatures > median(AUC_signatures),],
         aes(x=factor(sig,levels=sig[order(AUC_signatures)]), col=(sig == it_sig),
             y=AUC_signatures))+
    geom_bar(stat = "identity")+theme_bw()+
    theme(axis.text.x=element_text(angle = 45, hjust = 1))+labs(x='Signature', y='AUC exposure')+
    guides(col="none")
  # ggtitle(it_sig)+geom_label(aes(x=-0.5, y =0, label='...'), fill=NA, lty=NA)
  ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_AUC.pdf"),
         height = 2.5, width = 6)
  saveRDS(AUC_signatures, paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, "_AUC.RDS"))

}
# ggplot(current_exposures_normalised_melt[current_exposures_normalised_melt$Var1 == 1,],
#        aes(x=L1, y=value, group=Var1))+geom_line()+facet_wrap(.~Var2)
# 
# ggplot(current_exposures_normalised_melt[current_exposures_normalised_melt$Var1 == 2,],
#        aes(x=L1, y=value, group=Var1))+geom_line()+facet_wrap(.~Var2)
# 
# rowSums(current_exposures[[1]])
# rowSums(current_exposures[[2]])

##------------------------------------------------------------------------------------##
## looking at one sample, how does it change?
# current_exposures_melt$Var2 <- factor(current_exposures_melt$Var2, levels=levels_sigs)
# adj_mats <- lapply(1:max(current_exposures_melt$Var1), function(it_sample){
#   first_sample <- normalise_rw(as(dcast(current_exposures_melt[current_exposures_melt$Var1 == it_sample,] %>% dplyr::select(-Var1), L1~Var2, value.var = 'value')[,-1], 'matrix'))
#   image(first_sample)
#   
#   ## find graph of first transition
#   first_sample[1:2,]
#   
#   adj_mat <- matrix(0, ncol=length(vec_sigs), nrow = length(vec_sigs))
#   adj_mat[vec_sigs == it_sig, ] <- first_sample[2,]
#   adj_mat[, vec_sigs == it_sig ] <- first_sample[2,]
#   rownames(adj_mat) <- colnames(adj_mat) <- vec_sigs
#   
#   adj_mat
# })
# adj_mats_sum <- Reduce('+', adj_mats)
# 
# adj_mat[adj_mat>0] <- 1
# plot(igraph::graph_from_adjacency_matrix(adjmatrix = adj_mat), )
##------------------------------------------------------------------------------------##

##------------------------------------------------------------------------------------##

current_exposures[[1]]
current_exposures[[2]]

bleedingAUC <- lapply(vec_sigs, function(i) readRDS(paste0("../../../results/exploratory/bleeding_signatures/", i, add, "_AUC.RDS")))
names(bleedingAUC) <- vec_sigs
bleedingAUCrbind <- do.call('rbind', bleedingAUC)

pheatmap::pheatmap()
pheatmap::pheatmap(bleedingAUCrbind, cluster_cols = F, cluster_rows = F)

bleedingAUCrbind_nodiag <- bleedingAUCrbind
diag(bleedingAUCrbind_nodiag) <- NA

pheatmap::pheatmap(bleedingAUCrbind, cluster_cols = F, cluster_rows = F)

bleedingAUCrbindm <- melt(bleedingAUCrbind)
bleedingAUCrbindm2 <- melt(bleedingAUCrbind_nodiag)
bleedingAUCrbindm$grouped <- (rowSums(bleedingAUCrbind, na.rm = T)) > 0.6

ggplot(bleedingAUCrbindm, aes(x=value, col=Var1, label=Var1))+geom_density()+
  #geom_label(aes(x=max(value)))+
  guides(col='none')+theme_bw()+scale_x_continuous(trans = "log2")+
  facet_wrap(.~grouped)

rownames_as_col <- function(i) cbind.data.frame(sig=rownames(i), i)

cbind.data.frame(mean=rowSums(bleedingAUCrbind, na.rm = T),
                 sd=apply(bleedingAUCrbind, 1, sd, na.rm = T),
                 sig=rownames(bleedingAUCrbind))

bleedingAUCrbindmsum <- rownames_as_col(melt(apply(bleedingAUCrbind, 1, median, na.rm = T)))
bleedingAUCrbindmsum$sig <- factor(bleedingAUCrbindmsum$sig, levels=bleedingAUCrbindmsum$sig[order(bleedingAUCrbindmsum$value)])
ggplot(bleedingAUCrbindmsum,
       aes(x=sig, y=value))+
  geom_bar(stat = "identity")+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
  

head(bleedingAUCrbindm)
ggplot(cbind.data.frame(mean=rowSums(bleedingAUCrbind, na.rm = T),
                        sd=apply(bleedingAUCrbind, 1, sd, na.rm = T),
                        sig=rownames(bleedingAUCrbind)),
       aes(x=factor(sig, levels=sig[order(value)]), y=value))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes())


bleedingAUCrbindm$equality <- bleedingAUCrbindm$Var1 == bleedingAUCrbindm$Var2
ggplot(bleedingAUCrbindm, aes(x=factor(Var1, levels=levels(bleedingAUCrbindmsum$sig)),
                              y=value, col=equality))+geom_boxplot()+scale_y_continuous(trans = "log2")+
  theme_bw()+  theme(axis.text.x=element_text(angle = 45, hjust = 1))
bleedingAUCrbindm[bleedingAUCrbindm$equality,]

ggplot(bleedingAUCrbindm, aes(x=factor(Var1, levels=levels(bleedingAUCrbindmsum$sig)),
                              y=value, col=equality))+geom_boxplot()+#scale_y_continuous(trans = "log2")+
  theme_bw()+  theme(axis.text.x=element_text(angle = 45, hjust = 1))+labs(x='Signatures', y='AUC exposures')+
  guides(col='none')
ggsave(paste0("../../../results/exploratory/bleeding_signatures/all_signatures_AUC_boxplots.pdf"),
       height = 2.5, width = 10.5)
