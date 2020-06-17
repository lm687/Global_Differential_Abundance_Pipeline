#########################################################
################### Work in progress ####################
#########################################################

## Comparison of simulated results from the inferred parameters, for D-M vs simpler Multinomial
debug = FALSE
if(debug){
  rm(list = ls())
  setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
}
  
library(optparse)
library(rstan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(MCMCpack)
library(plyr)
library(CompSign)
# library(parallel)

source("3_analysis/helper/helper_analyse_posteriors.R")

vector_models = c('M', 'DM')
nits = 20000

donors = read.table("../data/restricted/pcawg/icgc-dataset-1591612699408/donor.tsv",
                    stringsAsFactors = FALSE, sep = "\t", header = TRUE)
files_donors = read.table("../data/restricted/pcawg/repository_1567600367.tsv",
                          stringsAsFactors = FALSE, sep = "\t", header = TRUE)

## Read in ROO objects
list_CT = suggested_list_ct()
#it_features = c('nucleotidesubstitution1', 'nucleotidesubstitution3', 'signatures')
it_features = c('signatures')
source("3_analysis/helper/load_ROO.R")

opt=list();
uuid_flder =  "../data/inference/"
opt$uuid_folder_DM = "../data/inference/"
opt$uuid_folder_M = "../data/inference/"
# opt$uuid_folder_LNM = "../data/inference/"

# fles = paste0(flder, opt$file_in)

# posteriorsDM =  "DM_PCAWG_Nits7000_signatures_PRAD-CA_f852b23b-fba9-4d39-8b72-e8c9dde15165.Rdata"
# posteriorsM = "M_PCAWG_Nits7000_signatures_PRAD-CA_bd605743-1434-408f-bd27-901c978ff354.Rdata"
# 
# posteriorsM = "M_PCAWG_Nits7000_signatures_LAML-KR_ccf08bf9-8675-4205-9345-aee2439456ef.Rdata"
# posteriorsDM =  "DM_PCAWG_Nits7000_signatures_LAML-KR_10afa802-2310-401c-b75b-a20871f54de4.Rdata"
# posteriorsDM =  "DM_PCAWG_Nits7000_signatures_LAML-KR_cd2d3c87-acca-4ef8-8c2b-f2f23184d62d.Rdata"


####################################################################################################
## Checking dimensions of posteriors (number of features should be the same in M and DM)
####################################################################################################
## check when this happens
all_fles = lapply(c(opt$uuid_folder_DM, opt$uuid_folder_M#, opt$uuid_folder_LNM
                    ), list.files)
# df = cbind(model=rep(c('DM', 'M', 'LNM'), sapply(all_fles, length)),
df = cbind(model=rep(vector_models, sapply(all_fles, length)),
           feature_type=sapply(unlist(all_fles), function(j) sapply(j, function(i) gsub("ROO.RData", "", strsplit(i, "_")[[1]][2]))) %>% unlist,
           ct=sapply(unlist(all_fles), function(j) sapply(j, function(i) strsplit(i, "_")[[1]][1])) %>% unlist)

t(table(df[df[,'feature_type'] == "signatures",'ct'], df[df[,'feature_type'] == "signatures", 'model']))

for(i in unique(df[,'ct'])){
  .x = (df[df[,'ct'] == i,c('model', 'feature_type')])
  print(i)
  print(try(table(.x[,1], .x[,2])))
}

which_both_signature = sapply(unique(df[,'ct']), function(i){
  subset_df = df[(df[,'feature_type'] == "signatures") & (df[,'ct'] == i),]
  if(!is.null(dim(subset_df))){
    all(vector_models %in% subset_df[,'model'])
  }else{
    FALSE
  }
})

plotcontour <- function(patient_idx, group_idx, model_name,  add=TRUE, true_contour, ...){
  if(group_idx == 1){
    patient_idx = patient_idx
  }else if(group_idx == 2){
    patient_idx = patient_idx +  npatientsx2/2
  }
  
  if(model_name == 'DM'){
    pcomp_idx = prcomp_all
  }else if(model_name == 'M'){
    pcomp_idx = prcomp_all_M
  }else if(model_name == 'LNM'){
    pcomp_idx = prcomp_all_LNM
  }else{
    stop('Wrong <model_name>')
  }
  groupsplit = splits_df[[model_name]]
  df_subset = (scale(groupsplit[[patient_idx]],
                     center = TRUE, scale = FALSE) %*% pcomp_idx$rotation)
  
  title = paste0(model_name, ' ', c('Early', 'Late')[group_idx])
  #paste0(c('Dirichlet-Multinomial', 'Multinomial')[model_name], ' ', c('Early', 'Late')[group_idx])
  if(true_contour){
    ## plot contour plot
    contour(kde2d(df_subset[,1],
                  df_subset[,2], n=50), drawlabels=FALSE, nlevels=10, add=add,  col='#5c5c5c', ...,
            main=title)
  }else{
    ## plot simulated points instead
    if(!add){
      plot(df_subset[,1],
           df_subset[,2], main=title,  pch=19, col=alpha('black', 0.02), ...)
    }else{
      points(df_subset[,1],
             df_subset[,2], main=title,  pch=19, col=alpha('black', 0.02),  ...)
    }
  }
  
}

plot_whole_contour = function(group_idx, model_name, true_contour=TRUE){
  if(model_name == 'DM'){
    projected_observed_idx = projected_observed
    lims = apply(rbind(prcomp_res, projected_observed_idx), 2, function(i) c(min(i, na.rm = TRUE), max(i, na.rm = TRUE)))
  }else if(model_name == "M"){
    projected_observed_idx = projected_observed_M
    lims = apply(rbind(prcomp_res_M, projected_observed_idx), 2, function(i) c(min(i, na.rm = TRUE), max(i, na.rm = TRUE)))
  }else if(model_name == 'LNM'){
    projected_observed_idx = projected_observed_LNM
    lims = apply(rbind(prcomp_res_LNM, projected_observed_idx), 2, function(i) c(min(i, na.rm = TRUE), max(i, na.rm = TRUE)))
  }else{
    stop('Wrong <model_name>')
  }
  
  if(length(projected_observed_idx) == 1 & is.na(projected_observed_idx)){
    plot(0, 0)
  }else{
    
    plotcontour(patient_idx = 1, group_idx = group_idx, model_name = model_name, add = FALSE,
                xlim = c(lims[,1]), ylim = c(lims[,2]), true_contour = true_contour)
    if(npatientsx2 > 2){
      for(patient_idx in 2:(npatientsx2/2)){
        plotcontour(patient_idx = patient_idx, group_idx = group_idx, model_name = model_name, add = TRUE, true_contour)
      }
    }
    
    for(patient_idx in 1:(npatientsx2/2)){
      if(group_idx == 1){
        points(projected_observed_idx[patient_idx,1], projected_observed_idx[patient_idx,2], col='red', pch=19)
      }else if(group_idx == 2){
        points(projected_observed_idx[patient_idx +  npatientsx2/2,1], projected_observed_idx[patient_idx +  npatientsx2/2,2],
               col='red', pch=19)
      }
    }
  }
}

# listfles = lapply(c(opt$uuid_folder_M, opt$uuid_folder_DM, opt$uuid_folder_LNM), function(uuid_flder){
listfles = lapply(vector_models, function(model){
  listfles = list.files(uuid_flder)
  if(length(listfles) == 0){
    listfles = list.files(uuid_flder) ## retry
  }
  if(model == 'M'){
    listfles = listfles[grepl(paste0(nits, 'ROO'), listfles)]
  }else if(model == 'DM'){
    listfles = listfles[grepl(paste0(nits, '_DMROO'), listfles)]
  }else if(model == 'DM'){
    listfles = listfles[grepl(paste0(nits, '_LNMROO'), listfles)]
  }
  
  return(listfles)
})

listfM = listfles[[1]]
listfDM = listfles[[2]]
# listfLNM = listfles[[3]]

####################################################################################################
## Simulate data with inferred parameters
####################################################################################################

#type_features = c('features1', 'signatures')
type_features="signatures"
ct_in_inference_results_list = list(unique(sort(as.character(sapply(list.files(opt$uuid_folder_DM), function(strng) strsplit(strng, '_')[[1]][1])))),
                                    unique(sort(as.character(sapply(list.files(opt$uuid_folder_M), function(strng) strsplit(strng, '_')[[1]][1])))))
names(ct_in_inference_results_list) = type_features

plot_bool = FALSE
for (type_feature in type_features){
  ct_in_inference_results = (ct_in_inference_results_list[[type_feature]])
  
  for(ct in ct_in_inference_results){
    # plots_return_0_2 = mclapply(ct_in_inference_results, function(ct){
    cat(ct, '\n')
    cat(type_feature, '\t', ct, '\n')
    
    posteriorsM = listfM[grepl(ct, listfM) & grepl(type_feature, listfM)][1]
    posteriorsDM = listfDM[grepl(ct, listfDM) & grepl(type_feature, listfDM)][1]
    # posteriorsLNM = listfLNM[grepl(ct, listfLNM) & grepl(type_feature, listfLNM)][1]
    
    # vector_models = c('DM', 'M', 'LNM')
    files_rdata = c(paste0(opt$uuid_folder_DM, posteriorsDM), paste0(opt$uuid_folder_M,posteriorsM)#, paste0(opt$uuid_folder_LNM, posteriorsLNM)
                    )
    posteriors = lapply(files_rdata,
                        function(f){
                          if(substr(f, nchar(f), nchar(f)) == "/" | basename(f) == "NA"){
                            ## no file
                            NA
                          }else{
                            print(f)
                            load(f)
                            x = tryCatch(extract(fit_stan))
                            if(is.null(x)){
                              ## no samples
                              NA
                            }else{
                              x
                            }
                          }
                        })
    
    
    ## Simulate with the total number of mutations for DM
    rowsums_toll = unlist(lapply(objects_sigs_per_CT[[type_feature]][[ct]], rowSums))
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
    
    ## Compare the coefficients beta
    ## since not all have been run for the same number of iterations, subset the posteriors
    lengths_beta = sapply(posteriors, function(i) if(length(i) == 1){if(is.na(i)){NA}} else{dim(i$beta)[3]})
    vector_models = vector_models[!is.na(lengths_beta)]
    posteriors = posteriors[!is.na(lengths_beta)]
    lengths_beta = lengths_beta[!is.na(lengths_beta)]
    if(length(lengths_beta) == 0){next}
    dim_beta = unique(lengths_beta); stopifnot(length(dim_beta) == 1)
    selected_rows = lapply(posteriors, function(i) sample(dim(i$beta)[1], 1000, replace = FALSE))
    posteriors_subset_beta = lapply(1:length(posteriors), function(idx_posterior) do.call('rbind', lapply(1:dim_beta, function(idx_feature) select_feature(df_with_slices = posteriors[[idx_posterior]]$beta, idx_select = idx_feature)[selected_rows[[idx_posterior]],])))
    sapply(posteriors_subset_beta, function(i) dim(i)[1])
    posteriors_subset_beta_intercept = do.call('cbind', lapply(posteriors_subset_beta, function(i) i[,1]))
    posteriors_subset_beta_slope = do.call('cbind', lapply(posteriors_subset_beta, function(i) i[,2]))
    colnames(posteriors_subset_beta_intercept) = colnames(posteriors_subset_beta_slope) = vector_models
    
    pdf(paste0('../results/comparison_models/beta_pairs_', ct, '_', type_feature, '.pdf'))
    if(dim(posteriors_subset_beta_intercept)[2] == 1){
      plot(0, 0, main='Only one model - no comparison')
    }else{
      par(mfrow=c(2,1))
      pairs(posteriors_subset_beta_intercept, main='Beta intercept pairs plot')
      pairs(posteriors_subset_beta_slope, main='Beta slope pairs plot')
    }
    dev.off()
    
    posteriors_subset_beta_slope
    dim(posteriors[[1]]$u) ## [nits, n, d-1]
    dim(posteriors[[1]]$beta) ### [nits, 2, d-1]
    dim(objects_sigs_per_CT[[type_feature]][[ct]][[1]]) ## [n,d]
    
    ## as many columns as patients
    dim(select_feature(posteriors[[1]]$u, 1))
    
    match_file = match(gsub("_active", "", rownames(objects_sigs_per_CT[[type_feature]][[ct]][[1]])),
                       sapply(files_donors$File.Name, function(i) strsplit(i, '[.]')[[1]][1]))
    rownames(objects_sigs_per_CT[[type_feature]][[ct]][[1]])
    files_donors[match_file,c('ICGC.Donor','File.Name') ]
    age_donors = donors[match(files_donors[match_file,]$ICGC.Donor, donors$icgc_donor_id),]
    age_donors = age_donors[!is.na(age_donors$icgc_donor_id),c('icgc_donor_id', 'donor_age_at_last_followup')]
    
    nfeatures = dim(posteriors[[1]]$u)[3]
    png(paste0("../results/link_clinical/age_u_", ct, '_', type_feature, '.png'),
        width = nfeatures*1.7, height = length(posteriors)*1.7, units = "in", res = 300)
    par(mfrow=c(length(posteriors),nfeatures))
    for(idx_model in 1:length(posteriors)){
      for(idx_feature in 1:nfeatures){
        if(length(rep(age_donors$donor_age_at_last_followup, each=dim(posteriors[[idx_model]]$u)[1])) != length(as.vector(select_feature(posteriors[[idx_model]]$u, idx_feature)))){stop()}
        plt_age = cbind(rep(age_donors$donor_age_at_last_followup, each=dim(posteriors[[idx_model]]$u)[1]),
                        as.vector(select_feature(posteriors[[idx_model]]$u, idx_feature)))
        plt_age = plt_age[!is.na(plt_age[,1]),]
        if(dim(plt_age)[1] > 0){
          plt_age = plt_age[sample(1:nrow(plt_age), 2000),]
          plot(plt_age)
        }else{
          plot(0,0)
        }
      }
    }
    dev.off()
    
    if(plot_bool){
      
      if(bool_data_avilable[1]){
        npatientsx2 = dim(posteriors[[1]]$alpha)[2]
        patient_idx = 1
        sim_counts = lapply(1:npatientsx2, function(patient_idx) normalise_cl(apply(t(apply(do.call('cbind', select_person(posteriors[[1]]$alpha, patient_idx)),
                                                                                            1, MCMCpack::rdirichlet, n=1)),
                                                                                    1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
        
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
        projected_observed = (scale(normalise_rw(do.call('rbind', objects_sigs_per_CT[[type_feature]][[ct]])),
                                    center = TRUE, scale = FALSE) %*% prcomp_all$rotation)[,1:2]
        
      }else{
        sim_counts = prcomp_all = prcomp_res = projected_observed = cols = NA
      }
      
      ## need to put the actual points
      ## now do the same for multinomial
      
      if(bool_data_avilable[2]){
        npatientsx2 = dim(posteriors[[2]]$theta)[2]
        
        ## Simulate with the total number of mutations for Multinomial
        sim_counts_M = lapply(1:npatientsx2, function(patient_idx) normalise_cl(apply(do.call('cbind', select_person(posteriors[[2]]$theta, patient_idx)),
                                                                                      1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
        
        sim_counts_M = do.call('rbind', sim_counts_M)
        sim_counts_M = sim_counts_M[! (colSums(apply(sim_counts_M, 1, is.na)) > 0),]
        
        cols_M = rep(1:npatientsx2, each=dim(sim_counts_M)[1]/npatientsx2)
        subset_M = unlist(lapply(unique(cols_M), function(i) sample(x = which(cols_M == i),
                                                                    size = 1000, #round(0.1*sum(cols_M == i)),
                                                                    replace = FALSE )))
        sim_counts_M = sim_counts_M[subset_M,]
        cols_M = cols_M[subset_M]
        prcomp_all_M = prcomp((sim_counts_M), scale. = FALSE, center=TRUE)
        prcomp_res_M = prcomp_all_M$x[,c(1,2)]
        projected_observed_M = (scale(normalise_rw(do.call('rbind', objects_sigs_per_CT[[type_feature]][[ct]])),
                                      center = TRUE, scale = FALSE) %*% prcomp_all_M$rotation)[,1:2]
        
      }else{
        sim_counts_M = prcomp_all_M = prcomp_res_M = projected_observed_M = cols_M = NA
      }
      
      
      # if(bool_data_avilable[3]){
      #   npatientsx2 = dim(posteriors[[3]]$theta)[2]
      #   
      #   ## Simulate with the total number of mutations for Logistic-normal multinomial
      #   sim_counts_LNM = lapply(1:npatientsx2, function(patient_idx) normalise_cl(apply(do.call('cbind', select_person(posteriors[[3]]$theta, patient_idx)),
      #                                                                                   1, rmultinom, n=1, size=rowsums_toll[patient_idx])))
      #   
      #   sim_counts_LNM = do.call('rbind', sim_counts_LNM)
      #   sim_counts_LNM = sim_counts_LNM[! (colSums(apply(sim_counts_LNM, 1, is.na)) > 0),]
      #   
      #   cols_LNM = rep(1:npatientsx2, each=dim(sim_counts_LNM)[1]/npatientsx2)
      #   subset_LNM = unlist(lapply(unique(cols_LNM), function(i) sample(x = which(cols_LNM == i),
      #                                                                   size = 1000,
      #                                                                   replace = FALSE )))
      #   sim_counts_LNM = sim_counts_LNM[subset_LNM,]
      #   cols_LNM = cols_LNM[subset_LNM]
      #   prcomp_all_LNM = prcomp((sim_counts_LNM), scale. = FALSE, center=TRUE)
      #   prcomp_res_LNM = prcomp_all_LNM$x[,c(1,2)]
      #   projected_observed_LNM = (scale(normalise_rw(do.call('rbind', objects_sigs_per_CT[[type_feature]][[ct]])),
      #                                   center = TRUE, scale = FALSE) %*% prcomp_all_LNM$rotation)[,1:2]
      #   
      # }else{
      #   sim_counts_LNM = prcomp_all_LNM = prcomp_res_LNM = projected_observed_LNM = cols_LNM = NA
      # }
      
      select_rows = function(df, colours){
        if(is.na(colours)){ NA}else{lapply(unique(colours), function(i) df[colours == i,])}
      }
      
      ## plotting the contours for all patients
      splits_df = list(select_rows(sim_counts, cols),
                       select_rows(sim_counts_M, cols_M)#,
                       #select_rows(sim_counts_LNM, cols_LNM)
                       )
      
      names(splits_df) = c('DM', 'M')#c('DM', 'M', 'LNM')
      # plot((scale(splits_df[[1]][[1]],
      #             center = TRUE, scale = FALSE) %*% prcomp_all$rotation)[,1:2], main='Dirichlet-Multinomial')
      # plot((scale(splits_df[[2]][[1]],
      #             center = TRUE, scale = FALSE) %*% prcomp_all_M$rotation)[,1:2], main='Multinomial')
      # 
      
      par(mfrow=c(1,3))
      # plot_whole_contour(group_idx = 1, model_name = 1, true_contour = FALSE)
      # plot_whole_contour(group_idx = 1, model_name = 1, true_contour = TRUE)
      
      png(paste0("../results/simulation_from_params/contourplots_", type_feature, "_", ct, ".png"))
      # par(mfrow=c(3,2))
      par(mfrow=c(2,2))
      plot_whole_contour(group_idx = 1, model_name = 'DM', true_contour = FALSE)
      plot_whole_contour(group_idx = 2, model_name = 'DM', true_contour = FALSE)
      plot_whole_contour(group_idx = 1, model_name = 'M', true_contour = FALSE)
      plot_whole_contour(group_idx = 2, model_name = 'M', true_contour = FALSE)
      # plot_whole_contour(group_idx = 1, model_name = 'LNM', true_contour = FALSE)
      # plot_whole_contour(group_idx = 2, model_name = 'LNM', true_contour = FALSE)
      dev.off()
    }
    # })
  }
}

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



type_features = c('features1', 'signatures')

ct_in_inference_results_list = list(as.character(sapply(list.files(opt$uuid_folder_DM), function(strng) strsplit(strng, '_')[[1]][5])),
                                    as.character(sapply(list.files(opt$uuid_folder_M), function(strng) strsplit(strng, '_')[[1]][5])))

names(ct_in_inference_results_list) = type_features
plot_bool = TRUE


fraction_in_credint = lapply(type_features, function(type_feature){
  ct_in_inference_results = (ct_in_inference_results_list[[type_feature]])
  
  fraction_in_credint = lapply(ct_in_inference_results, function(ct){
    cat(type_feature, '\t', ct, '\n')
    
    listfM = list.files(opt$uuid_folder_M)
    if(length(listfM) == 0){
      listfM = list.files(opt$uuid_folder_M) ## retry
    }
    listfDM = list.files(opt$uuid_folder_DM)
    if(length(listfDM) == 0){
      listfDM = list.files(opt$uuid_folder_DM) ## retry
    }
    posteriorsM = listfM[grepl(ct, listfM) & grepl(type_feature, listfM)]
    posteriorsDM = listfDM[grepl(ct, listfDM) & grepl(type_feature, listfDM)]
    
    files_rdata = c(paste0(opt$uuid_folder_DM, posteriorsDM), paste0(opt$uuid_folder_M,posteriorsM))
    posteriors = lapply(files_rdata,
                        function(f){
                          print(f)
                          load(f)
                          tryCatch(extract(fit_stan))
                        })
    
    ### Check overdispersion
    # posteriors[[1]] ## DM
    # posteriors[[2]] ## multinomial
    
    # comparison_overdispersion = mclapply(1:2, function(mdel_idx){
    
    nits = dim(posteriors[[2]]$theta)[1]
    size_subsample = 1000
    
    person_idx = 1
    group_idx = 1
    
    npeople = dim(posteriors[[2]]$theta)[2]/2 ## n people x 2
    stopifnot(npeople == nrow(objects_sigs_per_CT[[type_feature]][[ct]][[1]]))
    
    comparison_overdispersion = lapply(1:2, function(group_idx){
      x0 = lapply(1:npeople, function(person_idx){
        
        x1 = lapply(1:dim(posteriors[[2]]$theta)[3], function(signature_idx){
          cbind(rep(normalise_rw(objects_sigs_per_CT[[type_feature]][[ct]][[group_idx]])[person_idx,signature_idx], size_subsample),
                do.call('cbind', select_person(posteriors[[2]]$theta, person_idx+(group_idx-1)*npeople))[sample(1:nits, size_subsample,
                                                                                                                replace=FALSE),signature_idx])
        })
        do.call('rbind', x1)
      })
      do.call('rbind', x0)
    })
    
    length(unique(comparison_overdispersion[[2]][,1])) ## error
    
    comparison_overdispersion2 = do.call('rbind', comparison_overdispersion)
    comparison_overdispersion2[is.na(comparison_overdispersion2[,1]),1] = 0 ## these were zeros. turn into zeros
    # })
    
    png(paste0("../results/overdispersion/overdispersion_", paste0(ct, type_feature), ".png"))
    plot(comparison_overdispersion2[,1], comparison_overdispersion2[,2],
         xlab="True value", ylab = "Inferred value", main=paste0(ct, type_feature),
         col=rep(1:2, each=nrow(comparison_overdispersion2)/2), pch=19, cex=0.2)
    abline(0,1, lwd=3, lty=2, col='blue')
    dev.off()
    
    fctors = rep(1:(nrow(comparison_overdispersion2)/size_subsample), each=size_subsample)
    fraction_in_credint_bool = sapply(unique(fctors), function(idx){
      subst = comparison_overdispersion2[which(fctors == idx),]
      qntile = quantile(subst[,2], probs = c(0.25, 1-0.25))
      as.logical((subst[1,1] <= qntile[2]) &  (subst[1,1] >= qntile[1]))
    })
    
    return(sum(fraction_in_credint_bool)/length(fraction_in_credint_bool))
    
  })
  return(fraction_in_credint)
})

saveRDS(fraction_in_credint , file = "../data/robjects_cache/simulation_fraction_in_credint.RDS")

names(fraction_in_credint[[1]]) = ct_in_inference_results_list[[1]]
names(fraction_in_credint[[2]]) = ct_in_inference_results_list[[2]]
names(fraction_in_credint) = type_features

ggplot(melt(fraction_in_credint), aes(x=factor(L2, levels=names(sort(unlist(fraction_in_credint[[1]]), decreasing = TRUE))),
                                      y=value, shape=L1))+geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x="")
ggsave("../results/overdispersion/fraction_in_credint.pdf", width = 10)
fraction_in_credint = readRDS(".../data/robjects_cache/simulation_fraction_in_credint.RDS")

# apen = c()
# apen = c(apen, ("\\begin{figure}[h]\\centering"))
# for(ct in sort(unique(unlist(ct_in_inference_results_list)))){
#   if(ct %in% ct_in_inference_results_list[[1]]){
#     apen = c(apen, paste0("\\includegraphics[width=2in]{", "figures/overdispersion/overdispersion_", paste0(ct, 'features1'), ".png}\n"))
#   }
#   if(ct %in% ct_in_inference_results_list[[2]]){
#     apen = c(apen, (paste0("\\includegraphics[width=2in]{", "figures/overdispersion/overdispersion_", paste0(ct, 'signatures'), ".png}\n")))
#   }
# }
# apen = c(apen, ("\\end{figure}\n"))

cat(apen)
