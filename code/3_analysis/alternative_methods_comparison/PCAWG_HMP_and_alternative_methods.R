### Note: all of this is redundant. See Rmd instead.

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(Ternary)
library(MCMCpack)
library(TMB)
library(dplyr)
library(parallel)
library(umap)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

source("../../2_inference_TMB/helper_TMB.R")
source("../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

enough_samples = read.table("../../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
nonexogenous = read.table("../../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                          comment.char = "#", fill = F)
subset_sigs_sparse_cov_idx_nonexo <- read.table("../../../current/subset_sigs_sparse_cov_idx_nonexo.txt", stringsAsFactors = F, fill = T)

fles_roo <- list.files("../../../data/roo/", full.names = T)
fles_active <- fles_roo[grepl('_signaturesPCAWG_', fles_roo)]
fles_active <- fles_active[!grepl('subset_signatures', fles_active)]
signature_roo0 <- sapply(fles_active, readRDS)
names(signature_roo0) <- gsub("_signaturesPCAWG_ROO.RDS", "", basename(fles_active[grepl('_signaturesPCAWG_', fles_active)]))
signature_roo0 <- signature_roo0[match(enough_samples, names(signature_roo0))]
signature_roo_active <- lapply(signature_roo0, function(i) try(slot(i, 'count_matrices_active')))
names(signature_roo_active) <- enough_samples

# signature_roo0 = lapply(enough_samples, function(ct){
#   load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/", load_all_sigs = T)})
# signature_roo_active = lapply(enough_samples, function(ct){
#   load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../../data/")})

signature_roo_active_nonexo <- lapply(signature_roo_active, function(j){
  lapply(j, function(i) i[,!(colnames(i) %in% nonexogenous$V1)])
})

################################################################################################

diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMDL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/diagRE_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_M_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_M_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/fullRE_M_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_M <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/optim/diagRE_M_", ct, "_signatures.RDS")))
}, simplify = F)

sparseRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fulLRE_nonexo_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMDL_sortednonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fulLRE_sortednonexo_DM_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

fullRE_DMSL_SBS1 <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_SBS1baseline_DMSL_", ct, "signatures.RDS")))
}, simplify = F)

fullRE_halfDM <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_halfDM_", ct, "signatures.RDS")))
}, simplify = F)

fullRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/fullRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

diagRE_DMSL <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/diagRE_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

sparseRE_DMSL_nonexo <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_nonexo_DMSL_", ct, "_signatures.RDS")))
}, simplify = F)

################################################################################################

pvals_fullRE_M <- sapply(fullRE_M, function(i) try(wald_TMB_wrapper(i)))
pvals_fullRE_M <- p.adjust(pvals_fullRE_M)
pvals_diagRE_DM <- sapply(diagRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_diagRE_DM <- p.adjust(pvals_diagRE_DM)
pvals_DM <- sapply(sparseRE_DMSL, function(i) try(wald_TMB_wrapper(i)))
pvals_DM <- p.adjust(pvals_DM)
pvals_DMnonexo <- sapply(sparseRE_DMSL_nonexo, function(i) try(wald_TMB_wrapper(i)))
pvals_DMnonexo <- p.adjust(pvals_DMnonexo)
pvals_diagRE_DMSL_nonexo <- p.adjust(sapply(diagRE_DMSL_nonexo, function(i) try(wald_TMB_wrapper(i))))
################################################################################################


##2+ Sample DM-Distribution	
HMP_res_sigs <- lapply(signature_roo_active, function(signature_roo_it){
  try(HMP::Xdc.sevsample(list(t(signature_roo_it[[1]]),
                          t(signature_roo_it[[2]]))))
})

##2+ Sample Means w/o Reference Vector	
HMP_res_sigs_2 <- lapply(signature_roo_active, function(signature_roo_it){
  try(HMP::Xmcupo.sevsample(list(t(signature_roo_it[[1]]),
                              t(signature_roo_it[[2]]))))
})

HMP_res_estimates <- lapply(signature_roo_active, function(signature_roo_it){
  try(lapply(list((signature_roo_it[[1]]),
                              (signature_roo_it[[2]])), HMP::DM.MoM))
})
saveRDS(HMP_res_estimates, "../../../data/restricted/pcawg/pcawg_HMP_estimation_results.RDS")
HMP::Dirichlet.multinomial(Nrs = 30, shape = HMP_res_estimates$`Bone-Osteosarc`[[1]]$pi*HMP_res_estimates$`Bone-Osteosarc`[[1]]$theta)
length(HMP_res_estimates$`Bone-Osteosarc`[[1]]$pi)
(HMP_res_estimates$`Bone-Osteosarc`[[1]]$theta)
(HMP_res_estimates$`Bone-Osteosarc`[[2]]$theta)

HMP_res_sigs_nonexo <- lapply(signature_roo_active_nonexo, function(signature_roo_it){
  try(HMP::Xdc.sevsample(list(t(signature_roo_it[[1]]),
                              t(signature_roo_it[[2]]))))
})
names(HMP_res_sigs_nonexo) <- names(signature_roo_active)

plot(sapply(HMP_res_sigs, function(i) as.numeric(try(i$`p value`))),
     pvals_diagRE_DM)

plot(sapply(HMP_res_sigs_2, function(i) as.numeric(try(i$`p value`))),
     pvals_diagRE_DM)
## still doesn't make sense

plot(sapply(HMP_res_sigs, function(i) i$`p value`),
     log(pvals_diagRE_DM))

plot(sapply(HMP_res_sigs, function(i) i$`p value`),
     (pvals_fullRE_M))

plot(sapply(HMP_res_sigs_nonexo, function(i) as.numeric(try(i$`p value`))),
     pvals_DMnonexo)

plot(sapply(HMP_res_sigs_nonexo, function(i) as.numeric(try(i$`p value`))),
     pvals_diagRE_DMSL_nonexo)

###' It doesn't look similar at all
###' I need to show that, whenever the importance of random effects is high compared to fixed effects,
###' HMP is going to say that the datasets are not differentially abundant, when they are (there is
###' the case, e.g., of uterus adenocarcinoma

idx=23
names(HMP_res_sigs)[idx]
signature_roo_it= signature_roo_active[[idx]]
pvals_diagRE_DM[[idx]] ## differentially abundant
HMP_res_sigs[[idx]] ## not differentially abundant
grid.arrange(createBarplot(normalise_rw(signature_roo_it[[1]]))+guides(fill=FALSE),
createBarplot(normalise_rw(signature_roo_it[[2]]))+guides(fill=FALSE), nrow=1)

###' Can we get the coefficients from HMP?
HMP_res_sigs[[1]]
## doesn't look like it

