##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../../2_inference_TMB/helper_TMB.R")
source("../../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("../../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(mutSigExtractor)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
rel_path <- "../"
read_info <- function(ct){
  .x <- list(fullRE_M_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"))),
             # fullRE_DMSL_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"))),
             # fullRE_M_nonexo_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             # fullRE_DMSL_nonexo_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             # diagRE_DMDL_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             # diagRE_DMDL_nonexo_SP =  try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             # diagRE_DMDL_wSBS1SBS5nonexo_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMwSBS1SBS5nonexo_", ct, "_signaturesPCAWG.RDS"))),
             # fullREDMnoscaling_SP_nonexo =  try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", ct, "_signaturesPCAWG.RDS"))),
             # fullREDMnoscaling_SP_nonexo_subsets_and_amalgamations <- try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", ct, "_signaturesPCAWG.RDS"))),
             # fullREDMonefixedlambdanonexo_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             # fullREDMonefixedlambda2nonexo_SP = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambda2nonexo_", ct, "_signaturesPCAWG.RDS"))),
             # fullREDMonefixedlambdanonexo_SPSaA = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWGSaA.RDS"))),
             # fullREM_MSE = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesMSE.RDS"))),
             # fullREDM_MSE = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_signaturesMSE.RDS"))),
             # fullREDM_nucleotide1 = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDM_", ct, "_nucleotidesubstitution1.RDS"))),
             # diagREDM_MSE = try(readRDS(paste0(rel_path, "../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesMSE.RDS"))),
             # dataset_all_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = paste0(rel_path, "../../data/"), load_all_sigs = T),
             # dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = paste0(rel_path, "../../data/")),
             # dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = paste0(rel_path, "../../data/")),
             dataset_nucleotidesubstitution3 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3", path_to_data = paste0(rel_path, "../../data/"))
             # dataset_nucleotidesubstitution3MSE = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution3MSE", path_to_data =paste0(rel_path, "../../data/")),
             # dataset_active_sigs_MSE = load_PCAWG(ct = ct, typedata = "signaturesMSE", path_to_data = paste0(rel_path, "../../data/"), load_all_sigs = F),
             # DMM = list(z_DMM=lapply(1:8, function(k) try(read.table(paste0(rel_path, "../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.z"), sep = ',', skip = 1))),
             #            fit_DMM = lapply(1:8, function(k) try(read.table(paste0(rel_path, "../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.fit"), sep = ' '))))
  )
  # .x$dataset_nonexo <- give_subset_sigs_TMBobj(.x$dataset_active_sigs, nonexogenous$V1)
  # .x$dataset_nonexoSBS1SBS5 <- give_subset_sigs_TMBobj(.x$dataset_active_sigs, nonexogenouswSBS1SBS5$V1)
  # .x$colnames_notsorted_SP <- try(colnames(.x$dataset_active_sigs$Y))
  # .x$logR_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_notsorted_SP))
  # 
  # .x$colnames_nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexo$Y))
  # .x$logR_nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_nonexo_notsorted_SP))
  # 
  # ## nonexo, with SBS1 and SBS5
  # .x$colnames_wSBS1SBS5nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexoSBS1SBS5$Y))
  # .x$logR_wSBS1SBS5nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_wSBS1SBS5nonexo_notsorted_SP))
  return(.x)
}

read_info_list <- lapply(enough_samples, read_info)
names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##


## based on cosine similarity or other metrics, find which subsets of signatures fit the data best

created_files <- list.files("../../../results/exploratory/refitting_cosmic_sigs/")
created_files <- created_files[grepl('cossim_fall_plot.RDS', created_files)]

enough_samples <- enough_samples[ !(enough_samples %in% gsub("cossim_fall_plot.RDS", "", created_files)) ]

for(ct in enough_samples){
  out_RDS <- paste0("../../../results/exploratory/refitting_cosmic_sigs/", ct, "cossim_fall_plot.RDS")
  
  initial_nuc3 <- read_info_list[[ct]]$dataset_nucleotidesubstitution3
  initial_nuc3$Y
  
  remove_val <- function(rm_val, x){
    x[!(x == rm_val)]
  }
  current_sigs <- colnames(SBS_SIGNATURE_PROFILES_V3)
  
  cossims_mat <- matrix(NA, ncol(SBS_SIGNATURE_PROFILES_V3), ncol(SBS_SIGNATURE_PROFILES_V3)-1)
  rownames(cossims_mat) <- colnames(SBS_SIGNATURE_PROFILES_V3)
  best_removes <- rep(NA, ncol(SBS_SIGNATURE_PROFILES_V3) -1)
  
  for(it in 1:(ncol(SBS_SIGNATURE_PROFILES_V3) -1)){
    cat(it, '\n')
    current_fit <-  extract_sigs_TMB_obj(dataset_obj_trinucleotide=initial_nuc3,
                                         subset_signatures =  current_sigs)
    current_fit$Y
    cossiminit <- compare_signaturefit_to_data(tmb_obj_exposures = current_fit, initial_nuc3, SBS_SIGNATURE_PROFILES_V3, only_cosim=T)
    
    current_exposures <- colnames(current_fit$Y)
    
    all_possible_sig_removals <- lapply(current_exposures, function(sig_it){
      extract_sigs_TMB_obj(dataset_obj_trinucleotide=initial_nuc3,
                           subset_signatures = remove_val(sig_it, colnames(SBS_SIGNATURE_PROFILES_V3)))
    })
    
    all_possible_sig_removals
    
    cossims <- lapply(all_possible_sig_removals, function(j) compare_signaturefit_to_data(tmb_obj_exposures = j, initial_nuc3, SBS_SIGNATURE_PROFILES_V3, only_cosim=T))
      
    cossims_mat[current_exposures,it] <- unlist(cossims)
    
    best_remove <- which.min(cossiminit - unlist(cossims))
    best_removes[it] <- best_remove
    
    current_sigs <- current_sigs[-best_remove]
  }
  
  image(cossims_mat)
  
  max_cossims_mat_all <- apply(cossims_mat, 2, max, na.rm=T)
  plot(apply(cossims_mat, 2, max, na.rm=T))
  current_sigs
  
  "CNS-GBM"
  initial_nuc3$Y
  
  which_nas <- apply(cossims_mat, 2, function(i) which(is.na(i)))
  which_new_nas <- sapply(2:length(which_nas), function(j){
    if(j == 2){
      which_nas[[j]][!(which_nas[[j]] %in% colnames(SBS_SIGNATURE_PROFILES_V3))]
    }else{
      which_nas[[j]][!(which_nas[[j]] %in% which_nas[[j-1]])]
    }
    })
  
  which_new_nas <- c('All', names(which_new_nas), current_sigs)
  
  cossims_mat_losses <- apply(cossims_mat, 2, max, na.rm=T)
  cossims_mat_losses <- c(cossims_mat_losses, cossims_mat_losses[length(cossims_mat_losses)])
  ggplot(cbind.data.frame(cossim=cossims_mat_losses,
                          removed_sig = factor(which_new_nas, levels=which_new_nas)),
         aes(x=removed_sig, y=cossim, col=removed_sig == current_sigs))+geom_point()+geom_line(aes(group=1))+theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(x='Removed signature', y='Cosine similarity to observed trinucleotides')+
    ggtitle(ct)+guides(col=F)
  ggsave(paste0("../../../results/exploratory/refitting_cosmic_sigs/", ct, "cossim_fall_plot.pdf"), height = 3, width = 10)
  
  ## save object
  saveRDS(list(cossims_mat=cossims_mat, current_sigs=current_sigs, best_removes=best_removes, cossims_mat_losses=cossims_mat_losses,
               which_new_nas=which_new_nas), out_RDS)

}

new_refit_cossim <- extract_sigs_TMB_obj(dataset_obj_trinucleotide=read_info_list[[ct]]$dataset_nucleotidesubstitution3,
                     subset_signatures = c('SBS7b', 'SBS2', 'SBS8', 'SBS19', 'SBS1'))

give_barplot_from_obj(new_refit_cossim, legend_on = T)

colnames(read_info_list[[ct]]$dataset_active_sigs$Y)

# current_fitSBS1 <-  extract_sigs_TMB_obj(dataset_obj_trinucleotide=initial_nuc3,
#                                      subset_signatures =  current_sigs)
# 
# compare_signaturefit_to_data(tmb_obj_exposures = current_fitSBS1, initial_nuc3, SBS_SIGNATURE_PROFILES_V3, only_cosim=T)

# library(mutSigExtractor)
# sigs_strict <- fitToSignatures(
#   mut.context.counts=initial_nuc3$Y,
#   max.delta = 0.004,
#   detailed.output = F,
#   use.r.implementation = T,
#   signature.profiles=SBS_SIGNATURE_PROFILES_V3
# )

