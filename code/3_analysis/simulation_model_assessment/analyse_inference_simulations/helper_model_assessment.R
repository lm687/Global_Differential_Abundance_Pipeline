give_accuracies_with_varying_var <- function(var, two_var=F, datasets_arg=datasets, pvals_data_frame_arg=pvals_data_frame){
  if(two_var){
    do.call('rbind', apply(expand.grid(sapply(var, function(i) unique(unlist(sapply(datasets_arg, `[`, i))))), 1, function(vars_it){
      .pvals <- pvals_data_frame_arg[which(sapply(datasets_arg, '[', var[1]) == vars_it[[1]] & sapply(datasets_arg, '[', var[2]) == vars_it[[2]]),]
      .res_all_subset = put_vals_in_table(.pvals)
      .return <- cbind.data.frame(.res_all_subset, VAR1=vars_it[1], VAR2=vars_it[2], model=rownames(.res_all_subset))
      colnames(.return)[(ncol(.return)-2)] <- var[1]
      colnames(.return)[(ncol(.return)-1)] <- var[2]
      return(.return)
    }))    
  }else{
    do.call('rbind', lapply(unique(unlist(sapply(datasets_arg, `[`, var))), function(vars_it){
      .pvals <- pvals_data_frame_arg[which(sapply(datasets_arg, '[', var) == vars_it),]
      .res_all_subset = put_vals_in_table(.pvals)
      .return <- cbind.data.frame(.res_all_subset, d=vars_it, model=rownames(.res_all_subset))
      colnames(.return)[(ncol(.return)-1)] <- var
      return(.return)
    }))
  }
}

put_vals_in_table <- function(.pvals){
  warning('<fullREM> and <fullREDMSL> were interchanged until 20220201')
  rbind(fullREM=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_fullREM <= 0.05),
        fullREDMSL=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_fullREDMSL <= 0.05),
        diagREDMSL=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_diagREDMSL <= 0.05),
        diagREDM=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_diagREDM <= 0.05),
        pvals_chi_Harris=summarise_DA_detection(true = .pvals$true, predicted = .pvals$pvals_chi_Harris <= 0.05),
        ttest=summarise_DA_detection(true = .pvals$true, predicted = .pvals$ttest_props <= 0.05),
        ILR=summarise_DA_detection(true = .pvals$true, predicted = .pvals$ttest_ilr_adj <= 0.05),
        HMP=summarise_DA_detection(true = .pvals$true, predicted = .pvals$HMP <= 0.05),
        HMP2=summarise_DA_detection(true = .pvals$true, predicted = .pvals$HMP2 <= 0.05),
        perturbation=summarise_DA_detection(true = .pvals$true, predicted = .pvals$perturbation <= 0.05),
        permutation=summarise_DA_detection(true = .pvals$true, predicted = .pvals$permutation <= 0.05))
}

# colours_models <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
#                     fullREM= "#5B6DC8", HMP= "#3CA437", HMP2= "#6B244C" ,
#                     ILR= "#6ACDC5", permutation= "#DE1A1A" , perturbation= "#BBB53E",
#                     pvals_chi_Harris= "#2A297A", ttest=  "#995533"   )
# colours_models2 <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
#                     fullREM= "#5B6DC8", HMP= "#3CA437", HMP2= "#6B244C" ,
#                     ILR= "#6ACDC5", permutation= "#DE1A1A" , perturbation= "#BBB53E",
#                     chi_Harris= "#2A297A", ttest=  "#995533"   )
colours_models <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
                    fullREM= "#e69138", HMP= "#3CA437", HMP2= "#6B244C" ,
                    ILR= "#6ACDC5", permutation= "#8fce00" , perturbation= "#BBB53E",
                    pvals_chi_Harris= "#2A297A", ttest=  "#995533"   )
colours_models2 <- c(diagREDM= "#943CB4", diagREDMSL= "#194D44", fullREDMSL=  "#C6CF6E",
                     fullREM= "#e69138", HMP= "#3CA437", HMP2= "#6B244C" ,
                     ILR= "#6ACDC5", permutation= "#8fce00" , perturbation= "#BBB53E",
                     chi_Harris= "#2A297A", ttest=  "#995533"   )
