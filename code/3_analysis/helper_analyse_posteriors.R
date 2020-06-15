
plot_posterior_and_true = function(posterior, true, default_par=TRUE, include_other_lines=TRUE, ...){
  Nk = length(true)
  if(default_par){par(mfrow=c(1, length(true)))}
  
  if(Nk == 1){
    lims_param = t(rbind(min(c(posterior, true))-0.1,
                         max(c(posterior, true))+0.1))
    hist(posterior, xlim=c(lims_param[1],lims_param[2]), ...)
    abline(v=true, col='red')
  } else if (Nk > 1){
    sapply(1:Nk, function(k_it){
      lims_param = t(rbind(apply(rbind(posterior[,k_it], true[k_it]), 2, min)-0.1,
                           apply(rbind(posterior[,k_it], true[k_it]), 2, max)+0.1))
      hist(posterior[,k_it],
           # main=paste0('Posterior of part ', k_it), xlab="",
           # xlim=c(0,1),
           xlim=c(lims_param[k_it,1],lims_param[k_it,2]),
           ...);
      abline(v=true[k_it], col='red')
      if(include_other_lines)    for(i in 1:length(true)){abline(v=true[i], col='black')};
    })
  }
}

plot_posterior_and_true_generalisation = function(posterior, true, default_par=TRUE, ...){
  Nk = ncol(true)
  Ns = nrow(true)
  if(default_par){par(mfrow=c(Ns, Nk))}
  
  if(Nk == 1){
    stop('<plot_posterior_and_true_generalisation> with Nk == 1 needs to be written')
    # lims_param = t(rbind(min(c(posterior, true))-0.1,
    #                      max(c(posterior, true))+0.1))
    # hist(posterior, xlim=c(lims_param[1],lims_param[2]), ...)
    # abline(v=true, col='red')
  } else if (Nk > 1){
    sapply(1:Ns, function(n_it){
      sapply(1:Nk, function(k_it){
        lims_param = c(min(c(posterior[,,k_it][,n_it], true[n_it, k_it]))-0.1,
                       max(c(posterior[,,k_it][,n_it], true[n_it, k_it]))+0.1)
        hist(posterior[,,k_it][,n_it],
             # main=paste0('Posterior of part ', k_it), xlab="",
             # xlim=c(0,1),
             xlim=c(lims_param[1],lims_param[2]),
             ...);
        abline(v=true[n_it, k_it], col='red')
      })
    })
  }
}

get_bool_in_credible_interval = function(posterior, true){
  credible_int = quantile(posterior, probs = c(0.025, 0.975))
  (true >= credible_int[1]) && (true <= credible_int[2])
}

check_if_no_samples = function(posteriors, fles, flder){
  cat('Note! If they all appear as empty it might be because the fit as a different name (e.g. fit_DM, or fit_LNM)\nCheck this before removing the files.\n')
  sum(!(fles %in% names(posteriors)))
  ## remove empty files
  cat(paste0('for i in "', paste0(sapply(fles[! (fles %in% names(posteriors))], gsub, pattern = paste0(flder, '/', collapse = ''),
                                         replacement = "" ), collapse= '" "'), '"; do rm $i; done'))
}

posterior_covMat_to_Mean = function(posterior_cov_mat, d, mode=1){
  if(d == 2){
    if(mode == 1){
      mean(posterior_cov_mat[,,1])
    }else if(mode == 2){
      colMeans(posterior_cov_mat[,,1])
    }
  }else{
    cat('I am not sure if the matrix should be made row-wise or column-wise\n')
    sapply(1:(d-1), function(d_it)
      colMeans(posterior_cov_mat[,,d_it])
    )
  }
}


posterior_covMat_to_Median = function(posterior_cov_mat, d){
  if(d == 2){
    # median(posterior_cov_mat)
    apply(posterior_cov_mat[,,1], 2, median)
  }else{
    cat('I am not sure if the matrix should be made row-wise or column-wise\n')
    sapply(1:(d-1), function(d_it)
      apply(posterior_cov_mat[,,d_it], 2, median)
    )
  }
}


posterior_covMat_to_Mode = function(posterior_cov_mat, d){
  if(d == 2){
    # mode(posterior_cov_mat)
    apply(posterior_cov_mat[,,1], 2, mode)
  }else{
    cat('I am not sure if the matrix should be made row-wise or column-wise\n')
    sapply(1:(d-1), function(d_it)
      apply(posterior_cov_mat[,,d_it], 2, mode)
    )
  }
}

general_get_true_in_credint = function(fles, posteriors, true_params, name_param){
  mclapply(fles, function(f) outer(1:dim(posteriors[[f]][name_param][[1]])[3], 1:dim(posteriors[[f]][name_param][[1]])[2],
                                   Vectorize(function(idx_j, idx_l){
                                     get_bool_in_credible_interval(posterior = posteriors[[f]][name_param][[1]][,,idx_j][,idx_l],
                                                                   true = true_params[[f]][name_param][[1]][idx_l,idx_j])
                                   })))
}


modify_name = function(f, flder=flder, prefix="LNM_PCAWG") gsub(paste0(flder, prefix), "", f) %>% gsub(pattern = ".Rdata", replacement = "")%>%
  gsub(pattern = "_features1_", replacement = "")

shorten_name = function(names_str) sapply(names_str, function(name_str) paste0(strsplit(name_str, '-')[[1]][1:2], collapse ='-'))

## select person
select_person = function(df_with_slices, idx_select){
  lapply(1:dim(df_with_slices)[3], function(idx_slice) df_with_slices[,,idx_slice][,idx_select])
}

select_feature = function(df_with_slices, idx_select){
  df_with_slices[,,idx_select]
}

split_chains = function(posteriors, nchains, length_single_chain){
  dim_coefficient_isNA = F
  posteriors_list = list(); for(i in 1:nchains){posteriors_list[[i]] = posteriors}
  for(it_chains in 1:nchains){
    for(it_files in 1:length(posteriors)){
      for(it_features in names(posteriors_list[[it_chains]][[it_files]])){
        dim_coefficient = dim(posteriors_list[[it_chains]][[it_files]][[it_features]])[3]
        if(is.na(dim_coefficient)){dim_coefficient_isNA = T; dim_coefficient=1}
        for(it_slices in 1:dim_coefficient){
          # to_add = posteriors[[it_files]][[it_features]][,,it_slices][-c(1:additional_warmup),]
          interval_select_chains = ( 1+(length_chain*(it_chains-1))):(length_chain*it_chains)
          interval_select_chains = interval_select_chains[-c(1:additional_warmup)]
          if(dim_coefficient_isNA){
            posteriors_list[[it_chains]][[it_files]][[it_features]][-interval_select_chains] = NA
          }else{
            posteriors_list[[it_chains]][[it_files]][[it_features]][,,it_slices][-interval_select_chains,] = NA
          }
          dim_coefficient_isNA = F
        }
      }
      names(posteriors_list[[it_chains]][[it_files]]) = names(posteriors[[it_files]])
    }
    names(posteriors_list[[it_chains]]) = fles
  }
  names(posteriors_list) = paste0('Chain', 1:nchains)
  return(posteriors_list)
}

suggested_list_ct = function(){
#  sort(unique(gsub("_ROOSigs.RDS", "", list.files("/home/morril01/git_phd/out/ROO_PCAWG/deconvolution_clonal_subclonal_features1_ROO/"))))
  sort(unique(sapply(list.files("../data/roo/"), function(i) strsplit(i, '_')[[1]][1])))
}
