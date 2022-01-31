for(ct in enough_samples){
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
  ggsave(paste0("../../results/exploratory/refitting_cosmic_sigs/", ct, "cossim_fall_plot.pdf"), height = 3, width = 10)
  saveRDS(list(cossims_mat=cossims_mat, current_sigs=current_sigs, best_removes=best_removes, cossims_mat_losses=cossims_mat_losses,
               which_new_nas=which_new_nas),
          paste0("../../results/exploratory/refitting_cosmic_sigs/", ct, "cossim_fall_plot.RDS"))

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

