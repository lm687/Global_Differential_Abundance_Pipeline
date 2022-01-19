rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("recover_COSMIC_signatures.R")
source("../../2_inference_TMB/helper_TMB.R")
library(reshape2)
library(ggplot2)
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
    .x$exposuressigsQP <- as.numeric(.x$exposuressigsQP)
    .x <- dcast(.x[,c('exposuressigsQP', 'sig', 'sample')], sig~sample, value.var = 'exposuressigsQP')
    .x <- .x[,-1]
    current_exposures[[it]] <- t(.x)
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
  ggplot(current_exposures_melt,
         aes(x=L1, y=value, group=Var1))+geom_line(alpha=0.02)+
    facet_wrap(.~factor(Var2, levels=levels_sigs))+
    theme_bw()+labs(x='Iteration', y='Exposure')
  ggsave(paste0("../../../results/exploratory/bleeding_signatures/", it_sig, add, ".pdf"),
         height = 8, width = 8)
}
# ggplot(current_exposures_normalised_melt[current_exposures_normalised_melt$Var1 == 1,],
#        aes(x=L1, y=value, group=Var1))+geom_line()+facet_wrap(.~Var2)
# 
# ggplot(current_exposures_normalised_melt[current_exposures_normalised_melt$Var1 == 2,],
#        aes(x=L1, y=value, group=Var1))+geom_line()+facet_wrap(.~Var2)
# 
# rowSums(current_exposures[[1]])
# rowSums(current_exposures[[2]])

## looking at one sample, how does it change?
current_exposures_melt$Var2 <- factor(current_exposures_melt$Var2, levels=levels_sigs)
adj_mats <- lapply(1:max(current_exposures_melt$Var1), function(it_sample){
  first_sample <- normalise_rw(as(dcast(current_exposures_melt[current_exposures_melt$Var1 == it_sample,] %>% dplyr::select(-Var1), L1~Var2, value.var = 'value')[,-1], 'matrix'))
  image(first_sample)
  
  ## find graph of first transition
  first_sample[1:2,]
  
  adj_mat <- matrix(0, ncol=length(vec_sigs), nrow = length(vec_sigs))
  adj_mat[vec_sigs == it_sig, ] <- first_sample[2,]
  adj_mat[, vec_sigs == it_sig ] <- first_sample[2,]
  rownames(adj_mat) <- colnames(adj_mat) <- vec_sigs
  
  adj_mat
})
adj_mats_sum <- do.call(adj_mats, '+')

adj_mat[adj_mat>0] <- 1
plot(igraph::graph_from_adjacency_matrix(adjmatrix = adj_mat), )

