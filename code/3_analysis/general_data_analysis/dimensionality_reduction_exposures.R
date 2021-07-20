rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(umap)
library(ggplot2)
library(RColorBrewer)
source("../../2_inference_TMB/helper_TMB.R")

fles_roo <- list.files("../../../data/roo/", full.names = T)

enough_samples = readLines("~/Desktop/CT_sufficient_samples.txt")

signature_roo <- sapply(fles_roo[grepl('_signatures_', fles_roo)], readRDS)
signature_mutsigextractor_roo <- sapply(fles_roo[grepl('_signaturesmutSigExtractor_', fles_roo)], readRDS)

signature_roo <- sapply(signature_roo, function(i) try(slot(i, 'count_matrices_all')))
names(signature_roo) <- gsub("_signatures_ROO.RDS", "", basename(fles_roo[grepl('_signatures_', fles_roo)]))
signature_roo <- signature_roo[match(enough_samples, names(signature_roo))]
signature_roo <- signature_roo[sapply(signature_roo, typeof) == "list"]

signature_mutsigextractor_roo <- sapply(signature_mutsigextractor_roo, function(i) try(slot(i, 'count_matrices_all')))
names(signature_mutsigextractor_roo) <- gsub("_signaturesmutSigExtractor_ROO.RDS", "", basename(fles_roo[grepl('_signaturesmutSigExtractor_', fles_roo)]))
signature_mutsigextractor_roo <- signature_mutsigextractor_roo[match(enough_samples, names(signature_mutsigextractor_roo))]
signature_mutsigextractor_roo <- signature_mutsigextractor_roo[sapply(signature_mutsigextractor_roo, typeof) == "list"]


signature_roo_all <- do.call('rbind', lapply(signature_roo, function(i) rbind(i[[1]], i[[2]])))
signature_roo_all <- normalise_rw(signature_roo_all)
signature_roo_all_umap <- umap(signature_roo_all)

signature_roo_all_mutsigextractor <- do.call('rbind', lapply(signature_mutsigextractor_roo, function(i) rbind(i[[1]], i[[2]])))
signature_roo_all_mutsigextractor <- normalise_rw(signature_roo_all_mutsigextractor)
signature_roo_all_mutsigextractor_umap <- umap(signature_roo_all_mutsigextractor)


n <- length(signature_roo)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
names(col_vector) = names(signature_roo)

ggplot(cbind.data.frame(umap=signature_roo_all_umap$layout[,1:2], ct=rep(names(signature_roo), unlist(sapply(signature_roo, function(i) nrow(i[[1]])*2)))),
     aes(x=umap.1, y=umap.2, col=ct))+facet_wrap(.~ct)+
  geom_point(data=do.call('rbind', lapply(names(signature_roo), function(i) cbind.data.frame(umap=signature_roo_all_umap$layout[,1:2], ct=i))),
             col='gray', alpha=0.2)+
  geom_point()+
  scale_color_manual(values = col_vector)+guides(col=FALSE)+theme_bw()+labs(x='UMAP dimension #1',
                                                                            y='UMAP dimension #2')
ggsave("../../../results/exploratory/umap_exposures.pdf", height = 9, width = 9)


ggplot(cbind.data.frame(umap=signature_roo_all_mutsigextractor_umap$layout[,1:2], ct=rep(names(signature_mutsigextractor_roo), unlist(sapply(signature_mutsigextractor_roo, function(i) nrow(i[[1]])*2)))),
       aes(x=umap.1, y=umap.2, col=ct))+facet_wrap(.~ct)+
  geom_point(data=do.call('rbind', lapply(names(signature_mutsigextractor_roo), function(i) cbind.data.frame(umap=signature_roo_all_mutsigextractor_umap$layout[,1:2], ct=i))),
             col='gray', alpha=0.2)+
  geom_point()+
  scale_color_manual(values = col_vector)+guides(col=FALSE)+theme_bw()+labs(x='UMAP dimension #1',
                                                                            y='UMAP dimension #2')
ggsave("../../../results/exploratory/umap_exposures_mutsigextractor.pdf", height = 9, width = 9)



       