## change the name of files that didn't comnverge

setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/")
source("../2_inference_TMB/helper_TMB.R")

fles_roo <- list.files("../../data/roo/", full.names = T)

cts_for_PCAWG <- gsub("_signaturesPCAWG_ROO.RDS", "",
                      basename(fles_roo[grepl('_signaturesPCAWG_',
                                              fles_roo)]))
signatures_PCAWG <- sapply(cts_for_PCAWG,
                           load_PCAWG, typedata = "signaturesPCAWG",
                           path_to_data = "../../data/")
names(signatures_PCAWG) <- cts_for_PCAWG

length(signatures_PCAWG)
length(cts_for_PCAWG)

enough_samples <- unique(sapply(basename(fles_roo), function(i) strsplit(i, '_')[[1]][1]))

fullRE_M_SP <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS")))
}, simplify = F); names(fullRE_M_SP) <- enough_samples

for(ct in enough_samples[which(sapply(fullRE_M_SP, `[`, 'pdHess') == FALSE)]){
  system(paste0('mv ', paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"), ' ',
      paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct,  "_signaturesPCAWG_NC.RDS"), ''))
}

fullRE_DMSL_SP <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS")))
}, simplify = F); names(fullRE_DMSL_SP) <- enough_samples

for(ct in enough_samples[which(sapply(fullRE_DMSL_SP, `[`, 'pdHess') == FALSE)]){
  system(paste0('mv ', paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"), ' ',
                paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct,  "_signaturesPCAWG_NC.RDS"), ''))
}

fullRE_M_nonexo_SP <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS")))
}, simplify = F); names(fullRE_M_nonexo_SP) <- enough_samples

for(ct in enough_samples[which(sapply(fullRE_M_nonexo_SP, `[`, 'pdHess') == FALSE)]){
  system(paste0('mv ', paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"), ' ',
                paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct,  "_signaturesPCAWG_NC.RDS"), ''))
}



fullRE_DMSL_nonexo_SP <- sapply(enough_samples, function(ct){
  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS")))
}, simplify = F); names(fullRE_DMSL_nonexo_SP) <- enough_samples
for(ct in enough_samples[which(sapply(fullRE_DMSL_nonexo_SP, `[`, 'pdHess') == FALSE)]){
  system(paste0('mv ', paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"), ' ',
                paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct,  "_signaturesPCAWG_NC.RDS"), ''))
}
