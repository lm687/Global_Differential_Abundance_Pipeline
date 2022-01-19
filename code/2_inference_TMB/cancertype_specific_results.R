## Comments on the cancer-specific results of PCAWG

##-----------------------------------------------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
source("../2_inference_TMB/helper_TMB.R")
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")
source("../3_analysis/recovery_COSMIC_signatures/recover_COSMIC_signatures.R")

library(TMB)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(dplyr)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
enough_samples = read.table("../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]
enough_samples
nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t",
                          comment.char = "#", fill = F)
nonexogenouswSBS1SBS5 = read.table("../../data/cosmic/exogenous_signatures_SBS_withSBS1SBS5.txt", sep = "\t",
                                   comment.char = "#", fill = F)
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##
read_info <- function(ct){
  .x <- list(fullRE_M_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREM_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMSL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambda_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_M_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullRE_DMSL_nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMsinglelambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDM_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_nonexo_SP =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMnonexo_", ct, "_signaturesPCAWG.RDS"))),
             diagRE_DMDL_wSBS1SBS5nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/diagREDMwSBS1SBS5nonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMnoscaling_SP_nonexo =  try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMnoscaling_SP_nonexo_subsets_and_amalgamations <- try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMnoscalingnonexosubset_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambdanonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambda2nonexo_SP = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambda2nonexo_", ct, "_signaturesPCAWG.RDS"))),
             fullREDMonefixedlambdanonexo_SPSaA = try(readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/fullREDMonefixedlambdanonexo_", ct, "_signaturesPCAWGSaA.RDS"))),
             dataset_all_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/"),
             dataset_active_sigs = load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/"),
             dataset_nucleotidesubstitution1 = load_PCAWG(ct = ct, typedata = "nucleotidesubstitution1", path_to_data = "../../data/"),
             DMM = list(z_DMM=lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.z"), sep = ',', skip = 1))),
                        fit_DMM = lapply(1:8, function(k) try(read.table(paste0("../../data/roo_for_DMM_SPpcawg/DMM_output/", ct, "_signaturesPCAWG_all", k, "_dmm.fit"), sep = ' '))))
               )
  .x$dataset_nonexo <- give_subset_sigs_TMBobj(.x$dataset_all_sigs, nonexogenous$V1)
  .x$dataset_nonexoSBS1SBS5 <- give_subset_sigs_TMBobj(.x$dataset_all_sigs, nonexogenouswSBS1SBS5$V1)
  .x$colnames_notsorted_SP <- try(colnames(.x$dataset_all_sigs$Y))
  .x$logR_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_notsorted_SP))
  
  .x$colnames_nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexo$Y))
  .x$logR_nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_nonexo_notsorted_SP))
  
  ## nonexo, with SBS1 and SBS5
  .x$colnames_wSBS1SBS5nonexo_notsorted_SP <- try(colnames(.x$dataset_nonexoSBS1SBS5$Y))
  .x$logR_wSBS1SBS5nonexo_notsorted_SP <- try(vector_cats_to_logR(.x$colnames_wSBS1SBS5nonexo_notsorted_SP))
  return(.x)
}
##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

read_info_list <- lapply(enough_samples, read_info)
names(read_info_list) <- enough_samples
##-----------------------------------------------------------------------------------------------------##

plot_betas_scatter <- function(TMB_obj, names_cats=NULL, softmax=F){
  betas_mat <- t(matrix(python_like_select_name(TMB_obj$par.fixed, 'beta'), nrow=2))
  if(softmax){
    betas_mat <- t(softmax(t(rbind(betas_mat, 0))))
    names_cats <- c(sapply(names_cats, function(i) strsplit(i, '/')[[1]][1]),
                    strsplit(names_cats[1], '/')[[1]][2])
  }
  betas_df <- data.frame(betas_mat)
  if(!is.null(names_cats)){
    betas_df$logR = names_cats
  }else{
    betas_df$logR <- 1:nrow(betas_df)
  }
  a <- ggplot(betas_df, aes(x=X1, y=X2, label=logR))+labs(x='Beta intercept', y='Beta slope')
  if(!TMB_obj$pdHess){
    a <- a+geom_point(col='red')
  }else{
    a <- a+geom_point()
  }
  a <- a + geom_label_repel(size=3)+theme_bw()
  return(a)
}

give_mycolours <- function(n){
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  myColors <- col_vector[1:n]
}

plot_scatter <- function(m){
  m <- data.frame(m)
  names_sigs <- colnames(m)
  colnames(m) <- c('x1', 'x2')
  ggplot(m, aes(x=x1, y=x2))+geom_point()+theme_bw()+labs(x=names_sigs[1], y=names_sigs[2])
}

plot_covariance_mat <- function(TMB_obj, names_cats=NULL){
  .x <- L_to_cov(python_like_select_name(TMB_obj$par.fixed, 'cov_par_RE'),
                 d=length(python_like_select_name(TMB_obj$par.fixed, 'beta'))/2)
  if(!is.null(names_cats))  colnames(.x) <- rownames(.x) <- names_cats
  hm <- pheatmap::pheatmap(.x, main='Covariances')
  return(arrangeGrob(hm[[4]]))
}

give_plot_differential_precision <- function(TMB_obj){
  ovrdisp <- cbind.data.frame( plot_lambdas(TMB_obj, return_df=T, plot=F))
  if(TMB_obj$pdHess){
    col='red'
  }else{
    col='black'
  }
  
  ovrdisp[ovrdisp$name == 'Lambda 1','name'] = 'Early'
  ovrdisp[ovrdisp$name == 'Lambda 2','name'] = 'Late'
  
  ggplot(ovrdisp, aes(x=name,  y=`Estimate`))+
    geom_point(col=col)+
    geom_errorbar(aes(ymin=`Estimate`-`Std..Error`,
                      ymax=`Estimate`+`Std..Error`), width=.1)+
    theme_bw()+
    labs(x='', y='Log lambda')
}

give_DMM_plot <- function(TMB_obj){
  if(all(sapply(read_info_list[[ct]]$DMM$z_DMM, typeof) == 'character')){
    ggplot()+theme_bw()
  }else{
    ggplot(cbind(melt(normalise_rw(TMB_obj$dataset_all_sigs$Y)),
                 grouptiming=rep(c('Clonal','Subclonal'), each=dim(TMB_obj$dataset_all_sigs$Y)[1]/2),
                 groupDMM=paste0('DMM ', apply(TMB_obj$DMM$z_DMM[[which.min(unlist(sapply(TMB_obj$DMM$fit_DMM, `[`, 2)))]][,2:3], 1, which.max))),
           aes(x=Var1, y=value, fill=Var2))+geom_bar(stat = "identity")+facet_wrap(.~grouptiming+groupDMM, ncol=2)+
      scale_fill_manual(name = "grp",values = give_mycolours(ncol(read_info_list[[ct]]$dataset_all_sigs$Y)))+theme_bw()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }
}

logR_nucleotidesubs1 <- vector_cats_to_logR(sort(unlist( sapply(c('C', 'T'),
                                                                function(base){ paste0(base, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]) }))))

betas_df_to_softmax <- function(df_betas){
  .x <- softmax(c(df_betas$Estimate, 0))
  names(.x) <- c(gsub("/.*", "", df_betas$LogR), gsub("*./", "", df_betas$LogR[1]))
  .x
}

##-----------------------------------------------------------------------------------------------------##
## reports for all ct
for(ct in enough_samples){
  .a <- give_barplot_from_obj(read_info_list[[ct]]$dataset_all_sigs, legend_on = T, legend_bottom = T, plot=F, arg_title='')
  pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report1.pdf"), onefile=FALSE, height = 12)
  print(cowplot::plot_grid(cowplot::plot_grid(.a),
                     cowplot::plot_grid(plot_betas(read_info_list[[ct]]$fullRE_DMSL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = 'fullRE_DMSL_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = 'diagRE_DMDL_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_nonexo_SP, names_cats = read_info_list[[ct]]$logR_nonexo_notsorted_SP, title = 'diagRE_DMDL_nonexo_SP'),
                                        plot_betas(read_info_list[[ct]]$diagRE_DMDL_wSBS1SBS5nonexo_SP, names_cats = read_info_list[[ct]]$logR_wSBS1SBS5nonexo_notsorted_SP, title = 'diagRE_DMDL_wSBS1SBS5nonexo_SP')),
                     cowplot::plot_grid(plot_betas_scatter(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP),
                                        plot_betas_scatter(read_info_list[[ct]]$diagRE_DMDL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP, softmax=T)),
                     cowplot::plot_grid(cowplot::plot_grid(createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)),
                                                      order_labels = names(sort(rowSums(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)),
                                                                                decreasing = F)), remove_labels=T)+guides(fill=F)+ggtitle('Sorted by mut. toll'),
                                                      give_plot_differential_precision(read_info_list[[ct]]$diagRE_DMDL_SP), ncol=1),
                                        plot_covariance_mat(read_info_list[[ct]]$fullRE_DMSL_SP, names_cats = read_info_list[[ct]]$logR_notsorted_SP)),
                     nrow=4, rel_heights = c(2,2,1,2)))
  dev.off()
}
##-----------------------------------------------------------------------------------------------------##


##-----------------------------------------------------------------------------------------------------##

all_bleeding <- lapply(enough_samples, function(ct){
  .x <- give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                           abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                           rel_path = "../../", return_dataframe = T)
  .x$ct = ct
  .x
})

all_bleeding <- do.call('rbind', all_bleeding)

all_bleeding$bool = (all_bleeding$exposuressigsQP > all_bleeding$exposures)
all_bleeding$bool[all_bleeding$exposuressigsQP  == all_bleeding$exposures] <- NA ## neither increase nor decrease
table(is.na(all_bleeding$bool))
table(all_bleeding$bool)
# all_bleeding_summary$bool[all_bleeding_summary$bool] = 'Overestimate'
# all_bleeding_summary$bool[all_bleeding_summary$bool == 'FALSE'] = 'Underestimate'

all_bleeding_summary <- all_bleeding %>% group_by(sig, ct) %>% dplyr::summarise(mean_change=mean(bool, na.rm = T))
all_bleeding_summary$shape = all_bleeding_summary$mean_change > 0.65

library(jcolors)
library(viridis)
tikzDevice::tikz(file ="../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_signatures.tex",
                 width=7, height = 8)
ggplot((all_bleeding_summary), aes(x=ct, y=factor(sig, levels=rev(gtools::mixedsort(unique(all_bleeding_summary$sig)))),
                                   fill=mean_change))+geom_tile()+
  # scale_fill_viridis(begin = 0, end = 1)+
  scale_fill_gradientn(colours = c('red', 'yellow', 'blue'), values = c(0,0.5, 1))+
  # jcolors::scale_fill_jcolors_contin("pal3", reverse = TRUE, bias = 0.2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(data = all_bleeding_summary[all_bleeding_summary$shape == T,])+
  labs(x='Cancer types', y='Signatures', fill='Fraction of overestimates over changed exposures')+
  theme(legend.position = "bottom")+guides(fill=FALSE)
dev.off()

#ggsave(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/bleeding_signatures.pdf"), width=7)

##-----------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------##

# if(opt$model == "fullREM"){
# TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial.cpp",  "-std=gnu++17")
# dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_ME_multinomial"))
# # }else if(opt$model == "diagREM"){
# TMB::compile("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_ME_multinomial"))
# mod_model_name = "diagRE_M"
# # }else if(opt$model == "fullREDM"){
TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomial"))
# mod_model_name = "fullRE_DM"
# }else if(opt$model == "diagREDM"){
TMB::compile("../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
dyn.load(dynlib("../2_inference_TMB/mm_multinomial/diagRE_ME_dirichletmultinomial"))
# mod_model_name = "diagRE_DM"
# # }else if(opt$model =="fullREDMsinglelambda"){
# TMB::compile("../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("../2_inference_TMB/mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
# mod_model_name = "fullREDMsinglelambda"
# # use_nlminb=T
# # }else if(opt$model =="diagREDMsinglelambda"){
# TMB::compile("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
# mod_model_name = "diagREDMsinglelambda"
# # use_nlminb=T
# # }else if(opt$model =="FEDMsinglelambda"){
# TMB::compile("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/FE_dirichletmultinomial_single_lambda"))
# mod_model_name = "FEDMsinglelambda"
# # }else if(opt$model =="fullREDMnoscaling"){
# TMB::compile("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling.cpp",  "-std=gnu++17")
# dyn.load(dynlib("2_inference_TMB/mm_multinomial/fullRE_ME_dirichletmultinomialnoscaling"))
# mod_model_name = "fullREDMnoscaling"

fullRE_DMDL_features <- list()
amalgamated_extra <- list()
fullRE_DMDL_extra <- list()
diagRE_DMDL_extra <- list()
##-----------------------------------------------------------------------------------------------------##
## Bone-Osteosarc
ct <- "Bone-Osteosarc"

## in https://www.nature.com/articles/s41588-019-0390-2 they show two pops, one with and one without SBS3
createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)),
              order_labels = names(sort(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)[,'SBS3'],
                                        decreasing = F)), remove_labels=T)+ggtitle('Sorted by SBS3')
## run without the minor signatures


View(read_info_list[[ct]]$dataset_all_sigs$Y)
fullDM_no_small_sigs <- wrapper_run_TMB(object = give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_all_sigs,
                           sigs_to_remove = c('SBS13', 'SBS17a', 'SBS17b', 'SBS30')),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_no_small_sigs
wald_TMB_wrapper(fullDM_no_small_sigs, fail_non_converged = F)
fullDM_no_small_sigs_nonexo <- wrapper_run_TMB(object = give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_all_sigs,
                                                                         sigs_to_remove = unique(c(c('SBS13', 'SBS17a', 'SBS17b', 'SBS30'),
                                                                                            nonexogenous$V1))),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullDM_no_small_sigs_nonexo
wald_TMB_wrapper(fullDM_no_small_sigs_nonexo, fail_non_converged = F)

plot_scatter(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS17a', 'SBS17b')])
plot_scatter(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS2', 'SBS3')])


fullRE_DMDL_features =  wrapper_run_TMB(object = read_info_list[[ct]]$dataset_nucleotidesubstitution1,
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(cowplot::plot_grid(plot_betas(fullRE_DMDL_features, names_cats = c(vector_cats_to_logR(sort(unlist( sapply(c('C', 'T'),
                          function(base){ paste0(base, '>', c('A', 'C', 'G', 'T')[! (c('A', 'C', 'G', 'T') == base)]) })))))),
                          plot_betas(fullDM_no_small_sigs, names_cats = vector_cats_to_logR(colnames(give_subset_sigs_TMBobj(read_info_list[[ct]]$dataset_all_sigs,
                                                                                                                             sigs_to_remove = c('SBS13', 'SBS17a', 'SBS17b', 'SBS30'))$Y)), title = "Removing various sigs")),
                         give_DMM_plot(read_info_list[[ct]]),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"), nrow=3, rel_heights = c(1,2,2)))
dev.off()
##-----------------------------------------------------------------------------------------------------##
ct <- "Breast-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 15000)))


fullRE_DMDL_features[[ct]] =  wrapper_run_TMB(object = read_info_list[[ct]]$dataset_nucleotidesubstitution1,
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                        model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
fullRE_DMDL_nonhypermutated <- readRDS( file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))

saveRDS(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut), file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'ROO_nonhypermutated', '.RDS'))

reorder_cols_first_to_last_TMBobj <- function(TMB_obj){
  .x <- colnames(TMB_obj$Y)
  TMB_obj$Y <- cbind(TMB_obj$Y[,-1], TMB_obj$Y[,1])
  colnames(TMB_obj$Y) <- c(.x[-1], .x[1])
  TMB_obj
}

amalgamated_extra[[ct]] <- reorder_cols_first_to_last_TMBobj(give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                  list_groupings = c(list(c('SBS9', 'SBS17a', 'SBS40')),
                                                          as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS9', 'SBS17a', 'SBS40'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
diagRE_DMDL_extra[[ct]]

plot(density(read_info_list[[ct]]$dataset_all_sigs$Y[,'SBS37']))

amalgamated_extra[[paste0(ct, '_2')]] <- reorder_cols_first_to_last_TMBobj(give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS9', 'SBS17a', 'SBS40', 'SBS1')),
                                                                                                                  list(c('SBS5', 'SBS17b', 'SBS37', 'SBS41', 'SBS3')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS1', 'SBS9', 'SBS17a', 'SBS40', 'SBS5', 'SBS17b', 'SBS37', 'SBS41', 'SBS3'))]))))
diagRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)

plot_betas(diagRE_DMDL_extra[[paste0(ct, '_2')]], names_cats = vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)))

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(cowplot::plot_grid(plot_betas(fullRE_DMDL_features[[ct]],
                                                       names_cats = logR_nucleotidesubs1),
                                            plot_betas(fullRE_DMDL_nonhypermutated, names_cats = read_info_list[[ct]]$logR_notsorted_SP, title = "Removing hypermut samples"), rel_widths = c(2,3)),
                         cowplot::plot_grid(createBarplot(normalise_rw(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)),
                                                          order_labels = names(sort(non_duplicated_rows(read_info_list[[ct]]$dataset_all_sigs$Y)[,'SBS37'])),
                                                          remove_labels=T)+guides(fill=F)+ggtitle('Sorted by SBS37'),
                                            plot_betas(diagRE_DMDL_extra[[ct]],
                                    names_cats = vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='diagRE_DMDL amalgamation')),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         plot_grid(give_plot_bleeding(names_sigs =  c('SBS9', 'SBS40'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS9', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../", scales_x_free = T),
                         give_plot_bleeding(names_sigs =  c('SBS9', 'SBS17a', 'SBS40'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS9','SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../", scales_x_free = T), ncol=2),
                         nrow=4, rel_heights = c(1,2,2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-GBM"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 15000)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
saveRDS(object = diagRE_DMDL_nonhypermutated, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_nonhypermutated', '.RDS'))

plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP)

amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS11', 'SBS30')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS11', 'SBS30'))])))
fullRE_DMDL_nonhypermutated_amalgamated =  wrapper_run_TMB(object =give_subset_samples_TMBobj( amalgamated_extra[[ct]], samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(plot_betas(fullRE_DMDL_nonhypermutated_amalgamated, vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y))),
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS11'),
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                            rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS5', 'SBS30'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS5', 'SBS37'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS37'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"), nrow=1
                                            ),
                         nrow=3, rel_heights = c(1,2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-Medullo"

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
                         give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                                            abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                                            rel_path = "../../"),
                         cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS18'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"),
                                            give_plot_bleeding(names_sigs = c('SBS8', 'SBS30'),
                                                               abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                               rel_path = "../../"), nrow=1
                         ),
                         nrow=2, rel_heights = c(2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "CNS-PiloAstro"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 600)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                                                               list_groupings = c(list(c('SBS19', 'SBS23')),
                                                                                                                  as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS19', 'SBS23'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
wald_TMB_wrapper(fullRE_DMDL_extra[[ct]])

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y))),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS18'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS11'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS8', 'SBS30'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS30'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"), nrow=1
  # ),
  nrow=2, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "ColoRect-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 1e5)))
fullRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_nonhypermutated

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS28', 'SBS40', 'SBS45')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS28', 'SBS40', 'SBS45'))]))))
colSums(amalgamated_extra[[ct]]$Y)
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(fullRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                   abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                   rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS45', 'SBS40'),
                                      abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS45', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                                            rel_path = "../../"),
                   give_plot_bleeding(names_sigs = c('SBS28', 'SBS40'),
                                      abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS28', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                                                            rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Eso-AdenoCA"

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  # plot_betas(fullRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS40'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
  give_plot_bleeding(names_sigs = c('SBS28', 'SBS17b'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS28', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../")),
  nrow=3, rel_heights = c(2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Head-SCC"

strange_sample <- names(sort((read_info_list[[ct]]$dataset_all_sigs$Y[,'SBS33']), dec=T)[1])
diagRE_DMDL_nostrange_sample =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = strange_sample),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nostrange_sample

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  plot_betas(diagRE_DMDL_nostrange_sample, read_info_list[[ct]]$logR_notsorted_SP, title = "without c0812962-a345-48b7-aec0-01336c2d1eed"),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  
  plot_scatter(matrix(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_all_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_all_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS33'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS33'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS5', 'SBS16'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS16'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS3', 'SBS4'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS4'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS17a', 'SBS40'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../")),
  nrow=4, rel_heights = c(1, 2, 1, 1)))
dev.off()

plot_scatter(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS1', 'SBS33')])+scale_y_continuous(trans = "log2")+
  scale_x_continuous(trans = "log2")



##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-ChRCC"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 4000)))
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS17a', 'SBS17b', 'SBS5', 'SBS29')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS17a', 'SBS17b', 'SBS5', 'SBS29'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP, title='Without hypermutated sample'),
                     plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation without hypermutated sample'), nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # plot_scatter(matrix(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_all_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_all_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS17a'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS5', 'SBS17b'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS1', 'SBS17a'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS1', 'SBS17b'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-RCC.clearcell"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 10000)))
diagRE_DMDL_nonhypermutated =  wrapper_run_TMB(object = give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_nonhypermutated

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(diagRE_DMDL_nonhypermutated, read_info_list[[ct]]$logR_notsorted_SP, title='Without hypermutated sample'),
                     # plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation without hypermutated sample'),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # plot_scatter(matrix(read_info_list[[ct]]$dataset_all_sigs$Y[,c('SBS33')]/rowSums(read_info_list[[ct]]$dataset_all_sigs$Y), ncol=2)[!grepl(strange_sample, unique(rownames(read_info_list[[ct]]$dataset_all_sigs$Y))),])+labs(x='Clonal SBS33', y='Subclonal SBS33'),
  # cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS5', 'SBS17a'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS5', 'SBS17b'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS1', 'SBS17a'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17a'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../"),
  #                    give_plot_bleeding(names_sigs = c('SBS1', 'SBS17b'),
  #                                       abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS17b'), read_info_list[[ct]]$colnames_notsorted_SP)],
  #                                       rel_path = "../../")),
  nrow=3, rel_heights = c(1, 2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Kidney-RCC.papillary"

# kidney papillary needs SBS13 removed, because it's in extremely low quantity and probably many zeros

give_plot_bleeding_wrapper <- function(sigs_bleed, ...){
  give_plot_bleeding(names_sigs = sigs_bleed,
                   abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( sigs_bleed, read_info_list[[ct]]$colnames_notsorted_SP)],
                   rel_path = "../../", resort_sig_labels = F, ...)
}

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj( read_info_list[[ct]]$dataset_all_sigs,
                                                              list_groupings = c(list(c('SBS13', 'SBS2'),
                                                                                      c('SBS5', 'SBS29'),
                                                                                       c( 'SBS40', 'SBS22')),
    as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in%
                                                                  c('SBS13', 'SBS2', 'SBS22', 'SBS29', 'SBS5', 'SBS40'))]))))
colSums(amalgamated_extra[[ct]]$Y)
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_subset', '.RDS'))
wald_TMB_wrapper(fullRE_DMDL_extra[[ct]])

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS13', 'SBS5')),
                     give_plot_bleeding_wrapper( c('SBS13', 'SBS40')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS29', 'SBS5')),
                     give_plot_bleeding_wrapper( c('SBS29', 'SBS40')), ncol=2),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='Amalgamation'),
                     nrow=1),
  nrow=4, rel_heights = c(2, 1, 1, 2)))
dev.off()


##-----------------------------------------------------------------------------------------------------##

ct <- "Liver-HCC"  


colnames(dataset_subset$Y)
ncol(dataset_subset$Y)
colnames(read_info_list[[ct]]$dataset_all_sigs$Y)
ncol(read_info_list[[ct]]$dataset_all_sigs$Y)

typedata='signatures'
dataset_subset <- give_subset_sigs_TMBobj( read_info_list[[ct]]$dataset_all_sigs,
                  sigs_to_remove = c('SBS14', 'SBS17a', 'SBS17b', 'SBS18', 'SBS24', 'SBS28', 'SBS31', 'SBS35', 'SBS53', 'SBS54', 'SBS9'))
diagRE_DMDL_subset =  wrapper_run_TMB(object = dataset_subset,
                                               model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_subset
saveRDS(object = diagRE_DMDL_subset, file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_subset', '.RDS'))

## there are many active signatures
load_PCAWG(ct = ct, typedata = "signaturesPCAWG", path_to_data = "../../data/")
length(python_like_select_name(read_info_list[[ct]]$fullRE_M_SP$par.fixed, 'beta'))/2
length(python_like_select_name(diagRE_DMDL_subset$par.fixed, 'beta'))/2


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_subset, vector_cats_to_logR(colnames(dataset_subset$Y)), title='Subset of sigs'),
                     nrow=1),
  nrow=2, rel_heights = c(2, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Lung-SCC" 

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS1', 'SBS5', 'SBS8')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS1', 'SBS5', 'SBS8'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)))

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS1', 'SBS4'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS4'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"),
  give_plot_bleeding(names_sigs = c('SBS1', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS1', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"), ncol=2),
  nrow=2, rel_heights = c(2, 1)))
dev.off()

# par(mfrow=c(1,1))
# for(i in 1:nrow(sigs_cosmic)){
#   plot(as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[i,]), type='l', main=i)
# }
# 
# as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[39,]) ## high in SBS1 (but also quite high in another)
# as.numeric(cbind(sigs_cosmic[,c(1,5)], sigs_cosmic[,-c(1,5)])[35,]) ## high in SBS1 (but also quite high in another)
# 
# rownames(sigs_cosmic) <- gsub('>', '/', rownames(sigs_cosmic))
# 
# a <- read.table("../../data/restricted/pcawg/pcawg_restricted_snv/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/pcawg_restricted_snv_counts/0009b464-b376-4fbc-8a56-da538269a02f", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/0009b464-b376-4fbc-8a56-da538269a02f_cluster_assignments.txt.gz", sep='\t', header = T)
# a <- read.table("../../data/restricted/pcawg/", sep='\t', header = T)
# head(a)
# 
# a <- readRDS("../../data/roo/Biliary-AdenoCA_nucleotidesubstitution3_ROO.RDS")
# 
# plot(a[a$mutation == rownames(sigs_cosmic)[35],'VAF'])


##-----------------------------------------------------------------------------------------------------##

ct <- "Lymph-BNHL" 

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS13', 'SBS5', 'SBS56'),
                                                                                      c('SBS1', 'SBS40')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS13', 'SBS5', 'SBS56', 'SBS1', 'SBS40'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS56', 'SBS5'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS56', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS13', 'SBS5'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS13', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Subset of sigs'),
                     nrow=1),
  nrow=3, rel_heights = c(2, 1, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Lymph-CLL"


##-----------------------------------------------------------------------------------------------------##

ct <- "Ovary-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 10000)))

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS26'),
                                                                                      c('SBS41', 'SBS40')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS3', 'SBS26', 'SBS40', 'SBS41'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS26'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS26'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS3', 'SBS41'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS40', 'SBS26'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS40', 'SBS26'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"),
                     give_plot_bleeding(names_sigs = c('SBS40', 'SBS41'),
                                        abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS40', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                                        rel_path = "../../"), ncol=2),
  give_plot_bleeding(names_sigs = c('SBS5', 'SBS3', 'SBS41'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS5', 'SBS3', 'SBS41'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../"),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Amalgamation of sigs'),
                     nrow=1),
  nrow=5, rel_heights = c(1.5,0.5, 0.5, 0.5, 1)))
dev.off()


##-----------------------------------------------------------------------------------------------------##

ct <- "Panc-AdenoCA"
hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 21000)))

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS5', 'SBS17a', 'SBS20')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS3', 'SBS5', 'SBS17a', 'SBS20'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)

fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))
# fullRE_DMDL_extra[[ct]] = readRDS(paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))

amalgamated_extra[[paste0(ct, '_2')]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS3', 'SBS5')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS3', 'SBS5'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'diagRE_DMDL_nonhypermutated', '.RDS'))

amalgamated_extra[[paste0(ct, '_3')]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                                            list_groupings = c(list(c('SBS3', 'SBS5'),
                                                                                                    c('SBS8', 'SBS13')),
                                                                                               as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS3', 'SBS5', 'SBS8', 'SBS13'))]))))
diagRE_DMDL_extra[[paste0(ct, '_3')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_3')]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[paste0(ct, '_3')]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, '_3diagRE_DMDL_nonhypermutated', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)), title='Subset of sigs'),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS3', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS3', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS3', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS17a', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS17a', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS17a', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  rel_widths = c(1,1,1, 1), ncol=4),
  cowplot::plot_grid(give_plot_bleeding(names_sigs = c('SBS20', 'SBS5'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS20', 'SBS5'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
  give_plot_bleeding(names_sigs = c('SBS20', 'SBS40'),
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0))[match( c('SBS20', 'SBS40'), read_info_list[[ct]]$colnames_notsorted_SP)],
                     rel_path = "../../", resort=F),
                     give_plot_bleeding_wrapper( c('SBS51', 'SBS28')), ncol=3),
  nrow=4, rel_heights = c(0.5,.8,0.4, 0.4 )))
dev.off()

mean(normalise_rw(read_info_list[[ct]]$dataset_all_sigs$Y)[,'SBS3'])
mean(normalise_rw(read_info_list[["Panc-AdenoCA"]]$dataset_all_sigs$Y)[,'SBS3'])

mean(normalise_rw(read_info_list[[ct]]$dataset_all_sigs$Y)[,'SBS13'])
mean(normalise_rw(read_info_list[["Panc-AdenoCA"]]$dataset_all_sigs$Y)[,'SBS13'])


all_bleeding[(all_bleeding$ct == "Panc-AdenoCA") & (all_bleeding$sig == 'SBS3'),]

##-----------------------------------------------------------------------------------------------------##

ct <- "Panc-Endocrine"


pancreas <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                 PancAdenoCA = betas_df_to_softmax(plot_betas(read_info_list[["Panc-AdenoCA"]]$diagRE_DMDL_SP, read_info_list[["Panc-AdenoCA"]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')))

pancreas$PancAdenoCA <- pancreas$PancAdenoCA[match(names(pancreas$PancEndocrine), names((pancreas$PancAdenoCA)))]

pancreas2 <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                 PancAdenoCA = betas_df_to_softmax(plot_betas(diagRE_DMDL_extra[["Panc-AdenoCA"]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0("Panc-AdenoCA", '_2')]]$Y)), return_df = T) %>% dplyr::filter(type_beta == 'Slope')))
names(pancreas2$PancAdenoCA)[names(pancreas2$PancAdenoCA) == '3+'] <- '5'

pancreas2$PancAdenoCA <- pancreas2$PancAdenoCA[match(names(pancreas2$PancEndocrine), names((pancreas2$PancAdenoCA)))]

pancreas3 <- list(PancEndocrine=betas_df_to_softmax(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP, return_df = T) %>% dplyr::filter(type_beta == 'Slope')),
                  PancAdenoCA = betas_df_to_softmax(plot_betas(diagRE_DMDL_extra[[paste0("Panc-AdenoCA", "_3")]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0("Panc-AdenoCA", "_3")]]$Y)), return_df = T) %>% dplyr::filter(type_beta == 'Slope')))
names(pancreas3$PancAdenoCA)[names(pancreas3$PancAdenoCA) == '3+'] <- '5'
names(pancreas3$PancAdenoCA)[names(pancreas3$PancAdenoCA) == '8+'] <- '13'

pancreas3$PancAdenoCA <- pancreas3$PancAdenoCA[match(names(pancreas3$PancEndocrine), names((pancreas3$PancAdenoCA)))]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS5', 'SBS39')),
  give_plot_bleeding_wrapper( c('SBS5', 'SBS26')), ncol=2),
  cowplot::plot_grid(ggplot(cbind.data.frame(sig=names(pancreas$PancEndocrine),
                                             betaslopeest=do.call('cbind', pancreas),
                                             abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel()+theme_bw()+geom_smooth(method = lm)+
      scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F),
  ggplot(cbind.data.frame(sig=names(pancreas2$PancEndocrine), betaslopeest=do.call('cbind', pancreas2),
                          abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel(min.segment.length = 0.2, , max.overlaps = 100)+theme_bw()+geom_smooth(method = lm)+
    scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F),
  ggplot(cbind.data.frame(sig=names(pancreas3$PancEndocrine), betaslopeest=do.call('cbind', pancreas3),
                          abundance=(colSums(read_info_list[[ct]]$dataset_active_sigs$Y))),
         aes(x=betaslopeest.PancEndocrine, y=betaslopeest.PancAdenoCA, label=sig))+
    geom_point(aes(size=abundance))+geom_label_repel(min.segment.length = 0.2, max.overlaps = 100)+theme_bw()+geom_smooth(method = lm)+
    scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+guides(size=F), ncol=3),
  nrow=4, rel_heights = c(0.5,0.8, 0.5, 1)))
dev.off()


##-----------------------------------------------------------------------------------------------------##

ct <- "Prost-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 8000)))
hypermut

amalgamated_extra[[ct]] <- (give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS40', 'SBS33', 'SBS37')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS40', 'SBS33', 'SBS37'))]))))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  nrow=4, rel_heights = c(0.5,0.5, 1, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Skin-Melanoma.cutaneous"

amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj( read_info_list[[ct]]$dataset_all_sigs,
                                                              list_groupings = c(list(c('SBS5', 'SBS2', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d'),
                                                                                      c('SBS17a', 'SBS17b')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS5', 'SBS2', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS17a', 'SBS17b'))])))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
                     give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
                     give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  give_plot_bleeding_wrapper( c('SBS1', 'SBS7c')),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     nrow=1),
  nrow=7, rel_heights = c(0.5,1, .6, .6, 0.5, 0.5, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##
ct <- "Stomach-AdenoCA"

amalgamated_extra[[paste0(ct)]] <- (give_amalgamated_exposures_TMBobj(read_info_list[[ct]]$dataset_all_sigs,
                                                                            list_groupings = c(list(c('SBS20', 'SBS43', 'SBS44'),
                                                                                                    c('SBS26', 'SBS28', 'SBS9', 'SBS5')),
                                          as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS20', 'SBS43', 'SBS44', 'SBS26', 'SBS28', 'SBS5', 'SBS9'))]))))
diagRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct)]],
                                           model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
saveRDS(object = diagRE_DMDL_extra[[paste0(ct)]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, '_diagRE_DMDL_amalgamation', '.RDS'))


pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 12)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct)]]$Y)),
                                title='Sorted sigs diagRE_DMDL_SP amalgamation', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS3', 'SBS9', 'SBS20', 'SBS26', 'SBS28', 'SBS43', 'SBS44', 'SBS51')),
                     give_plot_bleeding_wrapper( c('SBS3', 'SBS9', 'SBS20', 'SBS26', 'SBS28', 'SBS43', 'SBS5', 'SBS51')),nrow=1),
  give_plot_bleeding_wrapper(c("SBS1",   "SBS2",   "SBS3",   "SBS5",   "SBS9",   "SBS13",  "SBS15",  "SBS17a", "SBS17b", "SBS18",
                               "SBS20",  "SBS21",  "SBS26",  "SBS28",  "SBS40" , "SBS41",  "SBS43",  "SBS44",  "SBS51",  "SBS58")),
  # cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
  #                               title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
  #                    nrow=1),
  nrow=5, rel_heights = c(0.5,0.5, 0.8, 0.5, 1)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Thy-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 3000)))
hypermut

amalgamated_extra[[ct]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                              list_groupings = c(list(c('SBS2', 'SBS13', 'SBS58')),
                                                                                 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS2', 'SBS13', 'SBS58'))])))
fullRE_DMDL_extra[[ct]] <- wrapper_run_TMB(object = amalgamated_extra[[ct]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[ct]]
saveRDS(object = fullRE_DMDL_extra[[ct]], file = paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/particular_runs/", ct, 'fullRE_DMDL_nonhypermutated', '.RDS'))



amalgamated_extra[[paste0(ct, '_2')]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                             list_groupings = c(list(c('SBS2', 'SBS13', 'SBS58', 'SBS40')),
                                                                                as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS2', 'SBS13', 'SBS58', 'SBS40'))])))
fullRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[paste0(ct, '_2')]]

c

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40', 'SBS58')),
  # give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a', 'SBS40')),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7a')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7b')), ncol=2),
  # cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7c')),
  #                    give_plot_bleeding_wrapper( c('SBS2', 'SBS5', 'SBS7d')), ncol=2),
  # give_plot_bleeding_wrapper( c('SBS1', 'SBS7c')),
  cowplot::plot_grid(plot_betas(fullRE_DMDL_extra[[ct]], vector_cats_to_logR(colnames(amalgamated_extra[[ct]]$Y)),
                                title='fullRE_DMDL Amalgamation of sigs', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_2')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)),
                                title='fullRE_DMDL Amalgamation 2 of sigs', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[paste0(ct, '_3')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_3')]]$Y)),
                                title='diagRE_DMDL Amalgamation 3 of sigs', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_3')]], vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_3')]]$Y)),
                                title='fullRE_DMDL Amalgamation 3 of sigs', sort_by_slope = T),
                     nrow=1),
  nrow=4, rel_heights = c(0.5,1, .6, .6)))
dev.off()

##-----------------------------------------------------------------------------------------------------##

ct <- "Uterus-AdenoCA"

hypermut <- unique(names(select_self(sort(rowSums(read_info_list[[ct]]$dataset_all_sigs$Y)) > 10000)))
hypermut

amalgamated_extra[[paste0(ct, '_nohypermut')]] <- give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut)

dim(read_info_list[[ct]]$dataset_all_sigs$Y)
dim(amalgamated_extra[["Uterus-AdenoCA_nohypermut"]]$Y)
diagRE_DMDL_extra[[paste0(ct)]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_nohypermut')]],
                                                         model = "diagRE_DM", use_nlminb=T, smart_init_vals=F)
diagRE_DMDL_extra[[paste0(ct)]]


amalgamated_extra[[paste0(ct, '_2')]] <- give_amalgamated_exposures_TMBobj(give_subset_samples_TMBobj( read_info_list[[ct]]$dataset_all_sigs, samples_to_remove = hypermut),
                                                                           list_groupings = c(list(c('SBS5', 'SBS28', 'SBS44', 'SBS26'),
                                                                                                   c('SBS40', 'SBS3', 'SBS6', 'SBS15', 'SBS14'),
                                                                                                   c('SBS10a', 'SBS10b'),
                                                                                                   c('SBS2', 'SBS13')),
 as.list(colnames(read_info_list[[ct]]$dataset_all_sigs$Y)[!(colnames(read_info_list[[ct]]$dataset_all_sigs$Y) %in% c('SBS5', 'SBS28', 'SBS44', 'SBS26', 'SBS6', 'SBS14', 'SBS15', 'SBS10a', 'SBS10b', 'SBS40', 'SBS3', 'SBS2', 'SBS13'))])))
colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)
fullRE_DMDL_extra[[paste0(ct, '_2')]] <- wrapper_run_TMB(object = amalgamated_extra[[paste0(ct, '_2')]],
                                           model = "fullRE_DM", use_nlminb=T, smart_init_vals=F)
fullRE_DMDL_extra[[paste0(ct, '_2')]]

pdf(paste0("../../results/results_TMB/pcawg/reports_per_cancer_type/", ct, "_report2.pdf"), onefile=FALSE, height = 15)
print(cowplot::plot_grid(
  cowplot::plot_grid(plot_betas(read_info_list[[ct]]$diagRE_DMDL_SP, read_info_list[[ct]]$logR_notsorted_SP,
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     nrow=1),
  cowplot::plot_grid(plot_betas(diagRE_DMDL_extra[[paste0(ct)]],
                                vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_nohypermut')]]$Y)),
                                title='Sorted sigs diagRE_DMDL_SP', sort_by_slope = T),
                     plot_betas(fullRE_DMDL_extra[[paste0(ct, '_2')]],
                                vector_cats_to_logR(colnames(amalgamated_extra[[paste0(ct, '_2')]]$Y)),
                                title='Amalgamation fullRE_DMDL_extra', sort_by_slope = T),
                                        nrow=1),
  give_plot_bleeding(names_sigs = read_info_list[[ct]]$colnames_notsorted_SP,
                     abundances = softmax(c(python_like_select_name(read_info_list[[ct]]$diagRE_DMDL_SP$par.fixed, 'beta')[c(T,F)], 0)),
                     rel_path = "../../"),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS28', 'SBS5', 'SBS40')),
                     give_plot_bleeding_wrapper( c('SBS3', 'SBS5', 'SBS40')), nrow=1),
  cowplot::plot_grid(give_plot_bleeding_wrapper( c('SBS6', 'SBS5', 'SBS40')),
                     give_plot_bleeding_wrapper( c('SBS15', 'SBS6')), nrow=1),
  nrow=5, rel_heights = c(0.5,0.5, 1, .5, .5)))
dev.off()

##-----------------------------------------------------------------------------------------------------##
source("../3_analysis/helper/pcawg.colour.palette.R")
pcawg_palette <- pcawg.colour.palette(x = gsub("\\..*", "", names(read_info_list)),  scheme = "tumour.subtype")
names(pcawg_palette) <- names(read_info_list)
sbs2_sbs8 <- lapply(read_info_list, function(i) try(normalise_rw(i$dataset_all_sigs$Y)[,c('SBS2', 'SBS8')]))
sbs13_sbs8 <- lapply(read_info_list, function(i) try(normalise_rw(i$dataset_all_sigs$Y)[,c('SBS13', 'SBS8')]))
head(melt(sapply(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"], function(i) i[,1]/i[,2])))
par(mfrow=c(1,1))
df_sbs2_sbs8 <- data.frame(do.call('rbind', sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"]))
df_sbs13_sbs8 <- data.frame(do.call('rbind', sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"]))
df_sbs2_sbs8$ct <- rep(names(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"]), sapply(sbs2_sbs8[sapply(sbs2_sbs8, typeof) == "double"], nrow))
df_sbs2_sbs8$clonal = unlist(sapply(sapply(table(df_sbs2_sbs8$ct), function(i) i/2), function(j) rep(c("Clonal", 'Subclonal'), each=j)))
df_sbs13_sbs8$ct <- rep(names(sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"]), sapply(sbs13_sbs8[sapply(sbs13_sbs8, typeof) == "double"], nrow))
df_sbs13_sbs8$clonal = unlist(sapply(sapply(table(df_sbs13_sbs8$ct), function(i) i/2), function(j) rep(c("Clonal", 'Subclonal'), each=j)))
cowplot::plot_grid(ggplot(df_sbs2_sbs8, aes(x=SBS2, y=SBS8, col=gsub("\\..*", "", ct)))+
                     geom_point(col='black', size=3, alpha=0.2)+
                     geom_point()+#scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
                     theme_bw()+
    facet_wrap(.~ct+clonal, nrow=2, scales = "free")+
    scale_color_manual(values = pcawg_palette)+
    geom_smooth(method = lm, col='black')+guides(col=F),
  ggplot(df_sbs13_sbs8, aes(x=SBS13, y=SBS8, col=gsub("\\..*", "", ct)))+
    geom_point(col='black', size=3, alpha=0.2)+
    geom_point()+#scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
    theme_bw()+
    facet_wrap(.~ct+clonal, nrow=2, scales = "free")+
  scale_color_manual(values = pcawg_palette)+
    geom_smooth(method = lm, col='black')+guides(col=F), nrow=2)
ggsave("../../results/results_TMB/pcawg/SBS8_APOBEC_correlation_raw.pdf", height = 5, width = 10)
##-----------------------------------------------------------------------------------------------------##
python_like_select(list.files("../../results/results_TMB/pcawg/reports_per_cancer_type"), ct)


