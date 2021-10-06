rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(compositions)
library(ggplot2)
library(ggrepel)
library(lsa)
library(pheatmap)

sigs_cosmic <- read.table(paste0("../../../data/cosmic/sigProfiler_SBS_signatures_2019_05_22.csv"),
                          stringsAsFactors = FALSE, sep = ',', header = TRUE)
sigs_cosmic
rownames(sigs_cosmic) <- paste0(sigs_cosmic$Type, '_', sigs_cosmic$SubType)
sigs_cosmic <- sigs_cosmic[,-c(1:2)]
colSums(sigs_cosmic)
sigs_cosmic <- sweep(sigs_cosmic, 2, colSums(sigs_cosmic), '/') ## they don't exactly sum to one: normalise
colSums(sigs_cosmic)

sigs_cosmic <- t(sigs_cosmic)

clr_sigs <- compositions::clr(sigs_cosmic)
pca <- prcomp(clr_sigs)
sum_pca <- summary(pca)

ggplot(cbind.data.frame(pca$x[,1:2], sig=rownames(clr_sigs)), aes(x=PC1, y=PC2, label=sig, col=sig))+geom_point()+geom_label_repel()+
  theme_bw()+
  labs(x=paste0('PC1 (', round(sum_pca$importance['Proportion of Variance',1]*100), '%)'),
       y=paste0('PC1 (', round(sum_pca$importance['Proportion of Variance',2]*100), '%)'))
ggsave("../../../results/exploratory/about_cosmic_sigs/pca_cosmic_sigs.pdf", width = 7.5, height = 4.8)

cosine_sims <- lsa::cosine(t(sigs_cosmic))

pdf("../../../results/exploratory/about_cosmic_sigs/cosinesims_cosmic_sigs.pdf", width =9.5, height = 9.5)
print(pheatmap::pheatmap(cosine_sims))
dev.off()

aitch_dist <- as(dist(clr_sigs), 'matrix')
pdf("../../../results/exploratory/about_cosmic_sigs/aitchdist_cosmic_sigs.pdf", width =9.5, height = 9.5)
print(pheatmap::pheatmap(aitch_dist))
dev.off()

pdf("../../../results/exploratory/about_cosmic_sigs/aitchdist_cosmic_sigs_noSBS23.pdf", width =9.5, height = 9.5)
print(pheatmap::pheatmap(aitch_dist[rownames(aitch_dist) != "SBS23", colnames(aitch_dist) != "SBS23"]))
dev.off()
