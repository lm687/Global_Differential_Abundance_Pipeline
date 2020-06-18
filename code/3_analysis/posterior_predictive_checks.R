#########################################################
################### Work in progress ####################
#########################################################

## Comparison of simulated results from the inferred parameters, for D-M vs simpler Multinomial

rm(list = ls())

library(optparse)
library(rstan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(MCMCpack)
library(plyr)
library(CompSign)
library(scales)
library(optparse)

debug = FALSE
if(debug){
  setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
  opt = list(); opt$files_posteriors = c("../data/inference/Kidney-RCC.papillary_signatures_20000_MROO.RData", "../data/inference/Kidney-RCC.papillary_signatures_20000_DMROO.RData")
}else{
  option_list = list(
    make_option(c("--files_posteriors"), type="character", default=NA, 
                help="File with the posterior, with directory included", metavar="character")
  );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  opt$files_posteriors = strsplit(opt$files_posteriors, " ")[[1]]
  print(opt$files_posteriors)
}

source("3_analysis/helper/helper_analyse_posteriors.R")
source("3_analysis/helper/helper_simulation.R")
files_posterior_split = sapply(opt$files_posteriors, function(i) strsplit(i, "_")[[1]])
print(files_posterior_split)
ct = unique(basename(files_posterior_split[1,]))
type_feature = unique(files_posterior_split[2,])
nits =  as.numeric(files_posterior_split[3,])
model = gsub("ROO.RData", "",files_posterior_split[4,])
names(nits) = model
donors = read.table("../data/restricted/pcawg/icgc-dataset-1591612699408/donor.tsv",
                    stringsAsFactors = FALSE, sep = "\t", header = TRUE)
files_donors = read.table("../data/restricted/pcawg/repository_1567600367.tsv",
                          stringsAsFactors = FALSE, sep = "\t", header = TRUE)

give_roo_wrapper = function(.it_features, .list_CT){
  it_features = .it_features
  list_CT = .list_CT
  source("3_analysis/helper/load_ROO.R", local = TRUE)
  return(objects_sigs_per_CT[[it_features]][[list_CT]])
}

ROO_object = give_roo_wrapper(type_feature, ct)


posteriors = lapply(opt$files_posteriors,
                    function(f){
                      if(substr(f, nchar(f), nchar(f)) == "/" | basename(f) == "NA" ){
                        # no file
                        NA
                      }else{
                        print(f)
                        load(f)
                        x = tryCatch(rstan::extract(fit_stan))
                        if(is.null(x)){
                          ## no samples
                          NA
                        }else{
                          x
                        }
                      }
                    })
names(posteriors) = model

## Simulate with the total number of mutations for DM
rowsums_toll = unlist(lapply(ROO_object, rowSums))
length(rowsums_toll) ## number of patients*2
  
bool_data_avilable = rep(TRUE, length(posteriors))
for(i in 1:length(posteriors)){
  if(is.null(posteriors[[i]])){
    ## either no samples in rstan object, or no file
    bool_data_avilable[i] = FALSE
  }else if(length(posteriors[[i]]) == 1){
    if(is.na(posteriors[[i]])){ bool_data_avilable[i] = FALSE }
  }
}
names(bool_data_avilable) = model
  
## since not all have been run for the same number of iterations, subset the posteriors
lengths_beta_all = sapply(posteriors, function(i) if(length(i) == 1){if(is.na(i)){NA}} else{dim(i$beta)[3]})

posteriors_all = posteriors
posteriors = posteriors[!is.na(lengths_beta_all)]
lengths_beta = lengths_beta_all[!is.na(lengths_beta_all)]
if(length(lengths_beta) == 0){ quit() } ## no posteriors to analyse
if(length(lengths_beta) == 1){bool_comparison=FALSE}else{bool_comparison=TRUE} ## if we are comparing models, we need at least twoz

####################################################################################################
##################################### Comparison of beta values ####################################
####################################################################################################
if(bool_comparison){
  dim_beta = unique(lengths_beta); stopifnot(length(dim_beta) == 1)
  selected_rows = lapply(posteriors, function(i) sample(dim(i$beta)[1], 1000, replace = FALSE))
  selected_rows_all = rep(NA, length(posteriors_all)); selected_rows_all[!is.na(lengths_beta_all)] = selected_rows
  posteriors_subset_beta = lapply(1:length(posteriors), function(idx_posterior) do.call('rbind', lapply(1:dim_beta, function(idx_feature) select_feature(df_with_slices = posteriors[[idx_posterior]]$beta, idx_select = idx_feature)[selected_rows[[idx_posterior]],])))
  posteriors_subset_beta_all = rep(NA, length(posteriors_all)); posteriors_subset_beta_all[!is.na(lengths_beta_all)] = posteriors_subset_beta
  posteriors_subset_beta_intercept = do.call('cbind', lapply(posteriors_subset_beta_all, function(i){if(is.na(i)){NA}else{i[,1]}}))
  posteriors_subset_beta_slope = do.call('cbind', lapply(posteriors_subset_beta_all, function(i){if(is.na(i)){NA}else{i[,2]}}))
  colnames(posteriors_subset_beta_intercept) = colnames(posteriors_subset_beta_slope) = model
  
  pdf(paste0('../results/comparison_models/beta_pairs_', ct, '_', type_feature, '.pdf'))
  if(sum(!(apply(posteriors_subset_beta_intercept, 2, function(i) all(is.na(i))))) == 1){
    plot(0, 0, main='Only one model - no comparison')
  }else{
    par(mfrow=c(2,1))
    pairs(posteriors_subset_beta_intercept, main='Beta intercept pairs plot')
    pairs(posteriors_subset_beta_slope, main='Beta slope pairs plot')
  }
  dev.off()
}

####################################################################################################
############################### Correlation of age with random effects #############################
####################################################################################################
dim(posteriors[[1]]$u) ## [nits, n, d-1]
dim(posteriors[[1]]$beta) ### [nits, 2, d-1]
dim(ROO_object[[1]]) ## [n,d]

match_file = match(gsub("_active", "", rownames(ROO_object[[1]])),
                   sapply(files_donors$File.Name, function(i) strsplit(i, '[.]')[[1]][1]))
rownames(ROO_object[[1]])
files_donors[match_file,c('ICGC.Donor','File.Name') ]
age_donors = donors[match(files_donors[match_file,]$ICGC.Donor, donors$icgc_donor_id),]
age_donors = age_donors[!is.na(age_donors$icgc_donor_id),c('icgc_donor_id', 'donor_age_at_last_followup')]

nfeatures = dim(posteriors[[1]]$u)[3]
png(paste0("../results/link_clinical/age_u_", ct, '_', type_feature, '.png'),
    width = 3*1.7, height = 3*length(posteriors)*1.7, units = "in", res = 300)
par(mfrow=c(length(posteriors),1))
for(idx_model in 1:length(posteriors)){
    if(length(rep(age_donors$donor_age_at_last_followup, each=dim(posteriors[[idx_model]]$u)[1])) != length(as.vector(posteriors[[idx_model]]$u))){stop()}
    plt_age = cbind(rep(age_donors$donor_age_at_last_followup, each=dim(posteriors[[idx_model]]$u)[1]),
                    as.vector(posteriors[[idx_model]]$u))
    plt_age = plt_age[!is.na(plt_age[,1]),]
    if(dim(plt_age)[1] > 0){
      plt_age = plt_age[sample(1:nrow(plt_age), 2000),]
      plot(plt_age)
    }else{
      plot(0,0)
    }
  }
dev.off()


#########################################################################################################
################# Lower-dimensional representation of posteriors and of observed values ################# 
#########################################################################################################
list_for_model = lapply(model, function(name_model){
  if(bool_data_avilable[name_model]){
    npatients = dim(posteriors_all[[name_model]]$u)[2]
    npatientsx2 = npatients*2
    patient_idx = 1
    if(name_model == 'DM'){
      sim_counts = lapply(1:npatientsx2, function(patient_idx) normalise_cl(apply(t(apply(do.call('cbind', select_person(posteriors_all[[name_model]]$alpha, patient_idx)),
                                                                                          1, MCMCpack::rdirichlet, n=1)),
                                                                                  1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
    }else{
      sim_counts = lapply(1:npatientsx2, function(patient_idx) normalise_cl(apply(do.call('cbind', select_person(posteriors[[name_model]]$theta, patient_idx)),
                                                                                      1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
    }
    
    sim_counts = do.call('rbind', sim_counts)
    
    sim_counts = sim_counts[! (colSums(apply(sim_counts, 1, is.na)) > 0),]
    
    par(mfrow=c(1,1))
    cols = rep(1:npatientsx2, each=dim(sim_counts)[1]/npatientsx2)
    subset = unlist(lapply(unique(cols), function(i) sample(x = which(cols == i),
                                                            size = 1000,#round(0.1*sum(cols == i)),
                                                            replace = FALSE )))
    sim_counts = sim_counts[subset,]
    cols = cols[subset]
    prcomp_all = prcomp(na.omit(sim_counts), scale. = FALSE, center=TRUE)
    prcomp_res = prcomp_all$x[,c(1,2)]
    projected_observed = (scale(normalise_rw(do.call('rbind', ROO_object)),
                                center = TRUE, scale = FALSE) %*% prcomp_all$rotation)[,1:2]
    
    }else{
      sim_counts = prcomp_all = prcomp_res = projected_observed = cols = npatientsx2 = NA
    }
    return(list(sim_counts=sim_counts, prcomp_all=prcomp_all, prcomp_res=prcomp_res, projected_observed=projected_observed, cols=cols, npatientsx2=npatientsx2))
  })
  
sim_counts = lapply(list_for_model, function(i) i$sim_counts)
prcomp_all = lapply(list_for_model, function(i) i$prcomp_all)
prcomp_all = lapply(list_for_model, function(i) i$prcomp_all)
prcomp_res = lapply(list_for_model, function(i) i$prcomp_res)
projected_observed = lapply(list_for_model, function(i) i$projected_observed)
cols = lapply(list_for_model, function(i) i$cols)
npatientsx2 = lapply(list_for_model, function(i) i$npatientsx2)
npatientsx2 = as.numeric(na.omit(unlist(unique(npatientsx2))))
if(length(npatientsx2) != 1){stop('The number of patients seems to be different across patients. Stopping.\n')}
npatients = npatientsx2/2
names(sim_counts) = names(prcomp_all) = names(prcomp_all) = names(prcomp_res) = names(projected_observed) = names(cols) = model

select_rows = function(df, colours){
  if(is.na(colours)){ NA}else{lapply(unique(colours), function(i) df[colours == i,])}
}

## plotting the contours for all patients
splits_df = lapply(1:length(sim_counts), function(i) select_rows(sim_counts[[i]], cols[[i]]) )

names(splits_df) = model

png(paste0("../results/simulation_from_params/contourplots_", type_feature, "_", ct, ".png"))
par(mfrow=c(length(model),2))
sapply(model, function(model_idx){
  plot_whole_contour(group_idx = 1, model_name = model_idx, true_contour = FALSE)
  plot_whole_contour(group_idx = 2, model_name = model_idx, true_contour = FALSE)
})
# plot_whole_contour(group_idx = 1, model_name = 'DM', true_contour = FALSE)
# plot_whole_contour(group_idx = 2, model_name = 'DM', true_contour = FALSE)
# plot_whole_contour(group_idx = 1, model_name = 'M', true_contour = FALSE)
# plot_whole_contour(group_idx = 2, model_name = 'M', true_contour = FALSE)
# plot_whole_contour(group_idx = 1, model_name = 'LNM', true_contour = FALSE)
# plot_whole_contour(group_idx = 2, model_name = 'LNM', true_contour = FALSE)
dev.off()

#########################################################################################################

#########################################################################################################
###### Plots showing whether the credible intervals of the posteriors capture the observed values  ######
#########################################################################################################
fraction_in_credint_bool = sapply(model, function(model_idx){
  size_subsample = 1e3
  comparison_overdispersion = lapply(1:2, function(group_idx){
    x0 = lapply(1:npatients, function(patient_idx){
      
      if(model_idx %in% c('M', 'LNM')){
        x1 = lapply(1:dim(posteriors[[model_idx]]$theta)[3], function(signature_idx){
          cbind(rep(normalise_rw(ROO_object[[group_idx]])[patient_idx,signature_idx], size_subsample),
                do.call('cbind', select_person(posteriors$M$theta, patient_idx+(group_idx-1)*npatients))[sample(1:nits[[model_idx]], size_subsample,
                                                                                                                replace=FALSE),signature_idx])
        })
      }else if(model_idx == 'DM'){
        x1 = lapply(1:dim(posteriors[[model_idx]]$alpha)[3], function(signature_idx){
          cbind(rep(normalise_rw(ROO_object[[group_idx]])[patient_idx,signature_idx], size_subsample),
                normalise_cl(apply(t(apply(do.call('cbind',
                                                   lapply(select_person(posteriors_all[[model_idx]]$alpha, patient_idx+(group_idx-1)*npatients),
                                                          function(i) i[sample(1:nits[[model_idx]], size_subsample, replace=FALSE)])),
                                                                                              1, MCMCpack::rdirichlet, n=1)),
                                                                                      1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
          
        
        })
      }
      do.call('rbind', x1)
    })
    unlist(x0) #do.call('rbind', x0)
  })
  
  
  comparison_overdispersion2 = do.call('rbind', comparison_overdispersion)
  comparison_overdispersion2[is.na(comparison_overdispersion2[,1]),1] = 0 ## these were zeros. turn into zeros
  
  fctors = rep(1:(nrow(comparison_overdispersion2)/size_subsample), each=size_subsample)
  fraction_in_credint_bool = sapply(unique(fctors), function(idx){
    subst = comparison_overdispersion2[which(fctors == idx),]
    qntile = quantile(subst[,2], probs = c(0.25, 1-0.25))
    as.logical((subst[1,1] <= qntile[2]) &  (subst[1,1] >= qntile[1]))
  })
  
  png(paste0("../results/overdispersion/overdispersion_", paste0(ct, type_feature), ".png"))
  plot(comparison_overdispersion2[,1], comparison_overdispersion2[,2],
       xlab="True value", ylab = "Inferred value", main=paste0(ct, type_feature),
       col=alpha(rep(1:2, each=nrow(comparison_overdispersion2)/2), 0.2), pch=19, cex=0.2)
  abline(0,1, lwd=3, lty=2, col='blue')
  dev.off()
  
  return(sum(fraction_in_credint_bool)/length(fraction_in_credint_bool))
})


saveRDS(paste0("../data/robjects_cache/fraction_in_credint_",gsub('.RData', '', ct) , ".RDS"))


# 
# type_features = c('features1', 'signatures')
# 
# ct_in_inference_results_list = list(as.character(sapply(list.files(uuid_flder), function(strng) strsplit(strng, '_')[[1]][5])),
#                                     as.character(sapply(list.files(uuid_flder), function(strng) strsplit(strng, '_')[[1]][5])))
# 
# names(ct_in_inference_results_list) = type_features
# plot_bool = TRUE
# 
# 
# fraction_in_credint = lapply(type_features, function(type_feature){
#   ct_in_inference_results = (ct_in_inference_results_list[[type_feature]])
# 
#   fraction_in_credint = lapply(ct_in_inference_results, function(ct){
#     cat(type_feature, '\t', ct, '\n')
# 
#     listfM = list.files(uuid_flder)
#     if(length(listfM) == 0){
#       listfM = list.files(uuid_flder) ## retry
#     }
#     listfDM = list.files(uuid_flder)
#     if(length(listfDM) == 0){
#       listfDM = list.files(uuid_flder) ## retry
#     }
#     posteriorsM = listfM[grepl(ct, listfM) & grepl(type_feature, listfM)]
#     posteriorsDM = listfDM[grepl(ct, listfDM) & grepl(type_feature, listfDM)]
# 
#     files_rdata = c(paste0(uuid_flder, posteriorsDM), paste0(uuid_flder,posteriorsM))
#     posteriors = posteriors_all
# 
#     ### Check overdispersion
#     # posteriors[[1]] ## DM
#     # posteriors[[2]] ## multinomial
# 
#     # comparison_overdispersion = mclapply(1:2, function(mdel_idx){
# 
#     nits = dim(posteriors[[2]]$theta)[1]
#     size_subsample = 1000
# 
#     person_idx = 1
#     group_idx = 1
# 
#     npeople = dim(posteriors[[2]]$theta)[2]/2 ## n people x 2
#     stopifnot(npeople == nrow(objects_sigs_per_CT[[type_feature]][[ct]][[1]]))
# 
#     comparison_overdispersion = lapply(1:2, function(group_idx){
#       x0 = lapply(1:npeople, function(person_idx){
# 
#         x1 = lapply(1:dim(posteriors[[2]]$theta)[3], function(signature_idx){
#           cbind(rep(normalise_rw(objects_sigs_per_CT[[type_feature]][[ct]][[group_idx]])[person_idx,signature_idx], size_subsample),
#                 do.call('cbind', select_person(posteriors[[2]]$theta, person_idx+(group_idx-1)*npeople))[sample(1:nits, size_subsample,
#                                                                                                                 replace=FALSE),signature_idx])
#         })
#         do.call('rbind', x1)
#       })
#       do.call('rbind', x0)
#     })
# 
#     length(unique(comparison_overdispersion[[2]][,1])) ## error
# 
#     comparison_overdispersion2 = do.call('rbind', comparison_overdispersion)
#     comparison_overdispersion2[is.na(comparison_overdispersion2[,1]),1] = 0 ## these were zeros. turn into zeros
#     # })
# 
#     png(paste0("../results/overdispersion/overdispersion_", paste0(ct, type_feature), ".png"))
#     plot(comparison_overdispersion2[,1], comparison_overdispersion2[,2],
#          xlab="True value", ylab = "Inferred value", main=paste0(ct, type_feature),
#          col=rep(1:2, each=nrow(comparison_overdispersion2)/2), pch=19, cex=0.2)
#     abline(0,1, lwd=3, lty=2, col='blue')
#     dev.off()
# 
#     fctors = rep(1:(nrow(comparison_overdispersion2)/size_subsample), each=size_subsample)
#     fraction_in_credint_bool = sapply(unique(fctors), function(idx){
#       subst = comparison_overdispersion2[which(fctors == idx),]
#       qntile = quantile(subst[,2], probs = c(0.25, 1-0.25))
#       as.logical((subst[1,1] <= qntile[2]) &  (subst[1,1] >= qntile[1]))
#     })
# 
#     return(sum(fraction_in_credint_bool)/length(fraction_in_credint_bool))
# 
#   })
#   return(fraction_in_credint)
# })
# 
# saveRDS(fraction_in_credint , file = "../data/robjects_cache/simulation_fraction_in_credint.RDS")
# 
# names(fraction_in_credint[[1]]) = ct_in_inference_results_list[[1]]
# names(fraction_in_credint[[2]]) = ct_in_inference_results_list[[2]]
# names(fraction_in_credint) = type_features
# 
# ggplot(melt(fraction_in_credint), aes(x=factor(L2, levels=names(sort(unlist(fraction_in_credint[[1]]), decreasing = TRUE))),
#                                       y=value, shape=L1))+geom_point()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="")
# ggsave("../results/overdispersion/fraction_in_credint.pdf", width = 10)
# fraction_in_credint = readRDS(".../data/robjects_cache/simulation_fraction_in_credint.RDS")
# 
# # apen = c()
# # apen = c(apen, ("\\begin{figure}[h]\\centering"))
# # for(ct in sort(unique(unlist(ct_in_inference_results_list)))){
# #   if(ct %in% ct_in_inference_results_list[[1]]){
# #     apen = c(apen, paste0("\\includegraphics[width=2in]{", "figures/overdispersion/overdispersion_", paste0(ct, 'features1'), ".png}\n"))
# #   }
# #   if(ct %in% ct_in_inference_results_list[[2]]){
# #     apen = c(apen, (paste0("\\includegraphics[width=2in]{", "figures/overdispersion/overdispersion_", paste0(ct, 'signatures'), ".png}\n")))
# #   }
# # }
# # apen = c(apen, ("\\end{figure}\n"))
# 
# cat(apen)


print_latex = T
if(print_latex){
  ## print for LaTeX
  ct_in_inference_results_list2 = ct_in_inference_results_list
  all_ct= sort(unique(unlist(ct_in_inference_results_list2)))
  for(ct in all_ct){
    featbool = FALSE; sigsbool = FALSE
    cat("\n\\begin{figure}[h]\n\\centering\n")
    if(ct %in% ct_in_inference_results_list2$features1){
      cat(paste0("\\includegraphics[width=.6\\textwidth]{figures/simulation_from_params/contourplots_features1_", ct, ".pdf}"), '\n')
      featbool = TRUE
    }
    if(ct %in% ct_in_inference_results_list2$signatures){
      cat(paste0("\\includegraphics[width=.6\\textwidth]{figures/simulation_from_params/contourplots_signatures_", ct, ".pdf}"), '\n')
      sigsbool = TRUE
    }
    if(featbool & sigsbool){cat(paste0('\\caption{Contour plots for ', ct, ' for features (above) and signatures (below), in PCA space.'))
    }else if(featbool){cat(paste0('\\caption{Contour plots for ', ct, ' for features, in PCA space.'))
    }else if(sigsbool){cat(paste0('\\caption{Contour plots for ', ct, ' for signatures, in PCA space.'))}
    cat(' The contour plots are drawn for simulated data with the inferred parameters. The observed data  are shown as red dots.}')
    cat("\\end{figure}")
    cat("\n\\newpage\n")
  }
}


