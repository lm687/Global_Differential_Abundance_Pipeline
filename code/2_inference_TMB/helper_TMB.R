softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

reflect_matrix = function(m){
  m[nrow(m):1,]
}

give_z_matrix = function(n_times_2){
  a = matrix(0, nrow = n_times_2/2, ncol = n_times_2/2)
  diag(a) = 1
  rbind(a, a)
}

give_x_matrix = function(n_times_2){
  cbind(rep(1, n_times_2), rep(c(0,1), each=n_times_2/2))
}

rsq = function (x, y) cor(x, y) ^ 2

load_PCAWG = function(ct, typedata, simulation=FALSE, path_to_data="../../data/", read_directly=FALSE, old_version_creating_X_Z=F){
  if(simulation | read_directly){
    fle = ct
    print(ct)
    cat('Reading file', fle, '\n')
  }else{
    fle = paste0(path_to_data, "roo/", ct, '_', typedata, "_ROO.RDS" )
  }
  objects_sigs_per_CT_features <- readRDS(fle)
  
  if(simulation){
    ##' get the first element of the list, because the second is the
    ##' set of the parameters used in the creation of the dataset
    objects_sigs_per_CT_features = objects_sigs_per_CT_features[[1]]
  }
  
  if(length(objects_sigs_per_CT_features) == 1){
    if(typeof(objects_sigs_per_CT_features) != "S4"){
      if(is.na(objects_sigs_per_CT_features)){
        return(NA)
      }
    }
  }
  
  if(typedata %in% c("nucleotidesubstitution1", "simulation")){
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
  }else if(typedata == "nucleotidesubstitution3"){
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) | (length(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) == 0)){
      ## no active signatures
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else{
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_active")
    }
    objects_sigs_per_CT_features = lapply(objects_sigs_per_CT_features, function(i){
      rwn = rownames(i)
      .x = apply(i, 2, as.numeric)
      rownames(.x) = rwn
      round(.x)
    })
  }else{
    stop('Check <typedata> argument')
  }
  
  d = ncol(objects_sigs_per_CT_features[[1]]) ## number of signatures or features
  n = nrow(objects_sigs_per_CT_features[[1]]) ## number of samples
  
  if(all(rownames(objects_sigs_per_CT_features[[1]]) == rownames(objects_sigs_per_CT_features[[2]]))){
    old_version_creating_X_Z = T
  }else{
    stop('Patients in input data are rearranged. Could not create matrices X and Z')
  }
  
  if(old_version_creating_X_Z){
    ## used directly up until 2 August 2021
    
    # matrix of fixed effects
    X = matrix(NA, nrow=2, ncol=2*n)
    X[1,] = 1
    X[2,] = rep(c(0,1), each=n)
    
    # covariate matrix
    Z0 = matrix(0, nrow=n, ncol=n)
    diag(Z0) = 1
    Z = t(rbind(Z0, Z0))
    
    ## The counts
    W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])
  }

  return(list(x=t(X), z=t(Z), Y=W))
}

# wrapper_run_TMB = function(ct, typedata, model, simulation=FALSE, allow_new_LNM=FALSE, object=NULL, sort_columns=F){
#   
#   if(!is.null(object)){
#     ## if the object of data and covariates is an argument
#     data = object
#   }else{
#     ## else, read from RDS file
#     if(!simulation){
#       cat(ct)
#       cat(typedata)
#     }
#     data = load_PCAWG(ct, typedata, simulation)
#     if(length(data) == 1){
#       if(is.na(data)){
#         return(warning('RDS object is NA'))
#       }
#     }
#     
#     if(sort_columns){
#       data$Y = data$Y[,order(colSums(data$Y), decreasing = T)]
#     }
#     
#   }
#   
#   data$Y = matrix(data$Y, nrow=nrow(data$Y))
#   data$x = (matrix(data$x, ncol=2))
#   
#   d <- ncol(data$Y)
#   n <- ncol(data$z) ## number of INDIVIDUALS, not samples
#   
#   if(model == "M"){
#     data$num_individuals = n
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_random_effects = matrix(rep(1, n)),
#       logSigma_RE=1
#     )
#     obj <- MakeADFun(data, parameters, DLL="ME_multinomial", random = "u_random_effects")
#   }else if(model == "fullRE_M"){
#     data$num_individuals = n
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial", random = "u_large")
#   }else if(model == "fullRE_Mcat"){
#     data$num_individuals = n
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial_categorical", random = "u_large")
#   }else if(model == "diagRE_M"){
#     data$num_individuals = n
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial", random = "u_large")
#   }else if(model == "DM"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_random_effects = matrix(rep(1, n)),
#       logSigma_RE=1,
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="ME_dirichletmultinomial", random = "u_random_effects")
#   }else if(model == "FE_DM"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="FE_dirichletmultinomial")
#   }else if(model == "fullRE_DM"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "sparseRE_DM"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_RE_part = 1,
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial_sparsecov", random = "u_large")
#     }else if(model == "fullRE_DMcat"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial_categorical", random = "u_large")
#   }else if(model == "fullRE_DM_singlelambda"){
#     data$num_individuals = n
# 
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#       log_lambda = 1
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_singlelambda_dirichletmultinomial", random = "u_large")
#   }else if(model == "diagRE_DM"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="diagRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "fullRE_DM_altpar"){
#     data$num_individuals = n
#     data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#     
#     parameters <- list(
#       beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                      nrow = 2, byrow=TRUE)),
#       u_large = matrix(rep(1, (d-1)*n), nrow=n),
#       logs_sd_RE=rep(1, d-1),
#       cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#       log_lambda = matrix(c(2,2))
#     )
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial_altpar", random = "u_large")
#   }else{
#     stop('Specify correct <model>\n')
#   }
#   obj$hessian <- TRUE
#   opt <- do.call("optim", obj)
#   opt
#   opt$hessian ## <-- FD hessian from optim
#   # obj$he()    ## <-- Analytical hessian
#   return(sdreport(obj))
# }

sort_columns_TMB = function(object){
  order_cats <- order(colSums(object$Y), decreasing = F) ## increasing so that the last category is the highest
  object$Y = object$Y[,order_cats]
  return(object)
}

sort_columns_TMB_SBS1 = function(object){
  if("SBS1" %in% colnames(object$Y)){
    order_cats <- c(which('SBS1' != colnames(object$Y)), which('SBS1' == colnames(object$Y)))
    object$Y = object$Y[,order_cats]
  }else{
    warning('There are no exposures for SBS1 in this sample. Keeping the same order')
  }
  return(object)
}


wrapper_run_TMB = function(model, object=NULL, smart_init_vals=T, use_nlminb=F, initial_params=NULL){
  ## sort_columns=F, 
  
  ## if the object of data and covariates is an argument
  data = object
  
  # if(sort_columns){
  #   order_cats = order(colSums(data$Y), decreasing = T)
  #   data$Y = data$Y[,order_cats]
  # }
  
  data$Y = matrix(data$Y, nrow=nrow(data$Y))
  data$x = (matrix(data$x, ncol=2))
  
  d <- ncol(data$Y) ## number of signatures
  n <- ncol(data$z) ## number of INDIVIDUALS, not samples
  
  if(smart_init_vals){
    require(nnet)
    .x_multinom = multinom(data$Y ~ data$x[,2])
    beta_init = t(coef(.x_multinom))
  }else{
    beta_init = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                        nrow = 2, byrow=TRUE))
  }
  
  parameters <- list(
    beta = beta_init,
    u_large = matrix(rep(1, (d-1)*n), nrow=n),
    logs_sd_RE=rep(1, d-1),
    cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
  )
  
  if(model == "fullRE_M"){
    data$num_individuals = n
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
    # )
    dll_name <- "fullRE_ME_multinomial"
    rdm_vec <- "u_large"
    # obj <- MakeADFun(data, parameters, DLL=, random = )
  }else if(model == "diagRE_M"){
    data$num_individuals = n
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
    # )
    dll_name <- "diagRE_ME_multinomial"
    rdm_vec <- "u_large"
    
    # obj <- MakeADFun(data, parameters, DLL=, random = "")
    }else if(model == "FE_DM"){
      data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))

      parameters$u_large = NULL
      parameters$logs_sd_RE = NULL
      parameters$cov_par_RE = NULL
      parameters <- list(parameters, log_lambda = matrix(c(2,2)))
      # parameters <- list(
      #   beta = beta_init,
      #   log_lambda = matrix(c(2,2))
      # )
      dll_name <- "FE_dirichletmultinomial"
      rdm_vec <- NULL
      
      # obj <- MakeADFun(data, parameters, DLL=)
  }else if(model == "fullRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- c(parameters,list(log_lambda = matrix(c(2,2))))
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
    #   log_lambda = matrix(c(2,2))
    # )
    dll_name <- "fullRE_ME_dirichletmultinomial"
    rdm_vec <- "u_large"
    # obj <- MakeADFun(data, parameters, DLL="", random = )
  }else if(model == "fullRE_DMonefixedlambda"){
    ## fixing one of the overdispersion parameters
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters,list(log_lambda = matrix(c(2))))
    dll_name <- "fullRE_ME_dirichletmultinomial_onefixedlambda"
    rdm_vec <- "u_large"
  }else if(model == "fullRE_DMonefixedlambda2"){
    ## fixing one of the overdispersion parameters
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters,list(log_lambda = matrix(c(2))))
    dll_name <- "fullRE_ME_dirichletmultinomial_onefixedlambda2"
    rdm_vec <- "u_large"
  }else if(model == "fullRE_DMonefixedlambda3"){
    ## fixing one of the overdispersion parameters
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters,list(log_lambda = matrix(c(2))))
    dll_name <- "fullRE_ME_dirichletmultinomial_onefixedlambda3"
    rdm_vec <- "u_large"
  }else if(model == "fullREDMnoscaling"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters, log_lambda = list(matrix(c(2,2))))
    parameters$logs_sd_RE = NULL
    print(parameters)
    dll_name <- "fullRE_ME_dirichletmultinomialnoscaling"
    rdm_vec <- "u_large"
  }else if(model == "diagRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters$log_lambda = matrix(c(2,2))
    parameters$cov_par_RE = NULL
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   log_lambda = matrix(c(2,2))
    # )
    dll_name <- "diagRE_ME_dirichletmultinomial"
    rdm_vec <- "u_large"
    # obj <- MakeADFun(data, parameters, DLL="", random = "")
  }else if(model  == "fullREDMsinglelambda"){
      data$num_individuals = n
      parameters$log_lambda = 1.1
      
      # parameters <- list(
      #   beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
      #                  nrow = 2, byrow=TRUE)),
      #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
      #   logs_sd_RE=rep(1, d-1),
      #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
      #   log_lambda = 1
      # )
      rdm_vec <- "u_large"
      dll_name <- "fullRE_dirichletmultinomial_single_lambda"
      # obj <- MakeADFun(data, parameters, DLL="", random = )
  }else if(model == "fullREhalfDM"){
  stop("fullREhalfDM used to be done with  fullRE_dirichletmultinomial_single_lambda")
  data$num_individuals = n
  
  parameters <- list(parameters, log_lambda = 1.1)
  rdm_vec <- "u_large"
  dll_name <- "CHANGETHIS"
  }else if(model == "fullREDMsinglelambda2"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- list(parameters, log_lambda = 1.1)
    
    # parameters <- list(
    #   beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
    #                  nrow = 2, byrow=TRUE)),
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
    #   log_lambda = 1
    # )
    rdm_vec <- "u_large"
    dll_name <- "fullRE_dirichletmultinomial_single_lambda2"
    # obj <- MakeADFun(data, parameters, DLL="", random = )
  }else if(model == "diagREDMsinglelambda"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- c(parameters, log_lambda = 1.1)
    parameters$cov_par_RE = NULL

    # Error below indicates   <parameters> list not being created correcty
    # Error in MakeADFun(data = Data, parameters = Parameters, random = Random)
    # :
    #   Only numeric matrices, vectors, arrays or factors can be interfaced
    # 
    # 
    
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   log_lambda = 2
    # )
    dll_name <- "diagRE_dirichletmultinomial_single_lambda"
    rdm_vec <- "u_large"
    
    # obj <- MakeADFun(data, parameters, DLL="", random = "")
  }else if(model == "FEDMsinglelambda"){
    data$num_individuals = NULL
    parameters$u_large = NULL
    parameters$logs_sd_RE = NULL
    parameters$cov_par_RE = NULL
    parameters <- list(parameters, log_lambda = 1.1)

    # parameters <- list(
    #   beta = beta_init,
    #   log_lambda = 2
    # )
    dll_name <- "FE_dirichletmultinomial_single_lambda"
    rdm_vec <- NULL
    
    # obj <- MakeADFun(data, parameters, DLL="")
  }else if(model == "fullRE_dirichletmultinomial_singlelambda_REv2"){
    data$num_individuals = n
    
    parameters$u_large = NULL
    parameters <- list(parameters,
                       u_large1 = matrix(runif(n)),
                       u_large2 = matrix(runif(n)),
                       u_large3 = matrix(runif(n)),
                       u_large4 = matrix(runif(n)),
                       log_lambda=1.2)
    
    # parameters <- list(
    #   beta = beta_init,
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
    #   log_lambda = 2,
    #   u_large1 = matrix(runif(n)),
    #   u_large2 = matrix(runif(n)),
    #   u_large3 = matrix(runif(n)),
    #   u_large4 = matrix(runif(n))
    # )
    rdm_vec <- c("u_large1", "u_large2","u_large3", "u_large4")
    dll_name <- "fullRE_dirichletmultinomial_single_lambda_REv2"
    # obj <- MakeADFun(data, parameters, DLL="",
    #                  random = )
  }else if(model == "fullRE_multinomial_REv2"){
    data$num_individuals = n
    parameters$u_large = NULL
    parameters <- list(parameters,
                       u_large1 = matrix(runif(n)),
                       u_large2 = matrix(runif(n)),
                       u_large3 = matrix(runif(n)),
                       u_large4 = matrix(runif(n)))
    # parameters <- list(
    #   beta = beta_init,
    #   logs_sd_RE=rep(1, d-1),
    #   cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
    #   u_large1 = matrix(runif(n)),
    #   u_large2 = matrix(runif(n)),
    #   u_large3 = matrix(runif(n)),
    #   u_large4 = matrix(runif(n))
    # )
    dll_name <- "fullRE_ME_multinomial_REv2"
    rdm_vec <- c("u_large1", "u_large2","u_large3", "u_large4")
    # obj <- MakeADFun(data, parameters, DLL="",
    #                  random = )
  }else{
    stop('Specify correct <model>\n')
  }
  
  if(!is.null(initial_params)){
    ## initial parameters are passed as arguments
    parameters <- initial_params
  }
  
  if(is.null(rdm_vec)){
    ## fixed effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name)
  }else{
    ## random effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name, random = rdm_vec)
  }
  
  if(use_nlminb){
    opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max)
  }else{
    obj$hessian <- TRUE
    opt <- do.call("optim", obj)
    opt
    opt$hessian ## <-- FD hessian from optim
  }
  return_report <- sdreport(obj)
  # if(sort_columns){
  #   ## return results in the original order
  #   order_cats
  #   return_report$par.fixed[grepl('beta', names(return_report$par.fixed))] = as.vector(matrix(python_like_select_name(return_report$par.fixed, 'beta'), nrow=2)[,order_cats])
  #   return_report$cov.fixed
  #   return_report
  # }
  
  return(return_report)
}

python_like_select = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = vector)]
}

python_like_select_name = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = names(vector))]
}

python_like_select_colnames = function(matrix, grep_substring){
  matrix[,grepl(pattern = grep_substring, x = colnames(matrix))]
}

python_like_select_rownames = function(matrix, grep_substring){
  matrix[grepl(pattern = grep_substring, x = rownames(matrix)),]
}

clean_name = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][2:3], collapse = ""))
}

clean_name_fullRE = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][3:4], collapse = ""))
}

clean_name_fullRE_2 = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][4:5], collapse = ""))
}

re_vector_to_matrix = function(vec_RE, dmin1){
  ## Random effects vector to matrix
  matrix(vec_RE, ncol=dmin1)
}

give_summary_per_sample = function(TMB_object, verbatim=T){
  if(is.null(TMB_object)){
    "Object doesn't exist"
  }else{
    if(typeof(TMB_object) %in% c("character", "logical")){
      return('Timeout or some error')
    }else{
      if(TMB_object$pdHess){
        return('Good')
      }else{
        return('Non-PD')
      }
    }
  }
}

give_summary_of_runs2 = function(vector_TMB_objects, long_return, verbatim=T){
  timeout_bool = sapply(vector_TMB_objects, typeof) == "character"
  hessian_positivedefinite_bool = sapply(vector_TMB_objects[!timeout_bool], function(i){
    if(length(i) == 1){FALSE}else{i$pdHess}})
  summary_runs = c(sum(timeout_bool), sum(!hessian_positivedefinite_bool), sum(hessian_positivedefinite_bool))
  names(summary_runs) = c('(failed) timeout', '(failed) non-positive semi-definite hessian', '(successful) positive semi-definite hessian')
  if(long_return){
    list(Timeout=names(vector_TMB_objects)[timeout_bool],
         hessian_nonpositivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][!hessian_positivedefinite_bool],
         hessian_positivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][hessian_positivedefinite_bool])
  }else{
    return(summary_runs)
  }
}


give_summary_of_runs = function(vector_TMB_objects, long_return, verbati=T){
  if(verbatim){
    stop('Check give_summary_of_runs2 instead')
  }
  timeout_bool = sapply(vector_TMB_objects, typeof) == "character"
  hessian_positivedefinite_bool = sapply(vector_TMB_objects[!timeout_bool], function(i) i$pdHess)
  summary_runs = c(sum(timeout_bool), sum(!hessian_positivedefinite_bool), sum(hessian_positivedefinite_bool))
  names(summary_runs) = c('(failed) timeout', '(failed) non-positive semi-definite hessian', '(successful) positive semi-definite hessian')
  if(long_return){
    list(Timeout=names(vector_TMB_objects)[timeout_bool],
         hessian_nonpositivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][!hessian_positivedefinite_bool],
         hessian_positivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][hessian_positivedefinite_bool])
  }else{
    return(summary_runs)
  }
}

get_summary_stan = function(model, typefeature){
  sapply(as.character(unique(samples_files$CT)), function(i){
    idx = which( (stan_results$CT == i) & (stan_results$features == typefeature) & (stan_results$model == model))
    if(length(idx) == 0){return("Object doesn't exist")}else{
      rw = stan_results[idx,]
      if(rw$divergent.transitions == "True"){
        return('Divergent transitions')
      } else if(rw$Cancelled..time.limit. == "True"){
        return('Timeout')
      } else if(!is.na(rw$Rhat.high) | !is.na(rw$ESS)){
        return('No good convergence')
      } else{
        return('Good')
      }
    }
  })
}

load_posteriors = function(fle_rdata){
  load(fle_rdata)
  list(posteriors_betas=posteriors_betas,
       posteriors_betas_slope=posteriors_betas_slope,
       num_not_containing_zero=num_not_containing_zero)
}

wald_generalised = function(v, sigma){
  # warning('20201218: sigma**(1/2) has now been replaced by (as we had before sometime in November) sigma')
  chisqrt_stat = t(v) %*% solve(sigma) %*% v
  pchisq(q = chisqrt_stat, df = length(v), lower.tail = FALSE)
}

wald_TMB_wrapper = function(i, verbatim=TRUE){
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = select_slope_2(which(names(i$par.fixed) == "beta"), verbatim=verbatim)
    if(!i$pdHess){
      ## didn't converge
      NA
    }else{
      if(length(idx_beta) == 1){
        if(is.na(idx_beta)){
          NA
        }else{
          if(verbatim) cat('Check data - slope appears to be of length one (binomial)\n')
          wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
        }
      }else{
        wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
      }
    }
  }
}


wald_TMB_wrapper_overdisp = function(i, verbatim=TRUE){
  ## wald test for the overdispersion parameter
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = which(names(i$par.fixed) == "log_lambda")
    if(!i$pdHess){
      ## didn't converge
      NA
    }else{
      scaled_coefs <- scale(i$par.fixed[idx_beta], center = T, scale = F)
      wald_generalised(v = as.vector(scaled_coefs),
                       sigma = i$cov.fixed[idx_beta,idx_beta])
    }
  }
}

ttest_TMB_wrapper_overdisp = function(i, verbatim=TRUE){
  ## wald test for the overdispersion parameter
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = which(names(i$par.fixed) == "log_lambda")
    if(!i$pdHess){
      ## didn't converge
      NA
    }else{
      .summary <- summary(i)
      loglambdas <- python_like_select_rownames(.summary, 'log_lambda')
      mean_coefs <- loglambdas[,1]
      SE_coefs <- loglambdas[,2]
      
      dmin1 <- nrow(python_like_select_rownames(.summary, 'beta'))/2
      num_patients <- nrow(python_like_select_rownames(.summary, 'u_large'))/dmin1
      tstatistic <- (mean_coefs[1]-mean_coefs[2])/sum(SE_coefs)
      df <- num_patients*2 - 2
      #' For significance testing, the degrees of freedom for this test
      #' is 2n − 2 where n is the number of participants in each group. 
      2*pt(-abs(tstatistic), df=df)
    }
  }
}


select_slope = function(i){
  stop('Deprecated due to error! use <select_slope_2> instead')
  i[(length(i)/2+1):(length(i))]
}

select_slope_2 = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    # if(verbatim) warning('As per 27 August 2020 it seems clear that this version, and not <select_slope>, is correct')
    i[c(F,T)]
  }else{
    i[,c(F,T)]
  }
}

select_intercept = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    i[c(T,F)]
  }else{
    i[,c(T,F)]
  }
}


vector_to_ct_list = function(vec){
  ##' given a vector which contains two or more types of features for the same cancer types,
  ##' convert it to a matrix with the pairing per cancer type
  which_sigs = (grep('signatures', names(vec)))
  which_nuc = (grep('nucleotidesubstitution1', names(vec)))
  if(sum(length(which_sigs)+length(which_nuc)) < length(vec)){stop('There is a third feature category')}
  ct_naked1 = gsub('signatures', '', names(vec)[which_sigs])
  ct_naked2 = gsub('nucleotidesubstitution1', '', names(vec)[which_nuc])
  ct_naked1[match(ct_naked1, ct_naked2)]
  ret = cbind(vec[which_sigs][match(ct_naked1, ct_naked2)], vec[which_nuc])
  rownames(ret) = gsub('signatures', '', rownames(ret))
  colnames(ret) = c('signatures', 'nucleotidesubsitution1')
  return(ret)
}

give_UNSTRUCTURED_CORR_t_matrix = function(vec, dim_mat){
  # #https://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html
  m = matrix(1, nrow = dim_mat, ncol = dim_mat)
  ## fill in the order that TMB's UNSTRUCTURED_CORR_t saves the covariances
  m[unlist(sapply(2:nrow(m), function(rw) seq(from = rw,length.out = (rw-1), by = nrow(m) )))] = vec
  m[unlist(sapply(2:nrow(m), function(cl) seq(from = (((cl-1)*nrow(m))+1),length.out = (cl-1), by = 1 )))] = vec
  return(m)
}

get_count_object = function(ct, feature, pre_path=NULL){
  if(is.null(pre_path)){
    pre_path = "../../data/roo/"
  }
  readRDS(paste0(pre_path, ct, "_", feature, "_ROO.RDS"))
}

get_count_object_file = function(fle){
  readRDS(fle)
}

give_stderr = function(i, only_slopes=T, only_betas=T){
  if(!only_betas){
    stop('Not yet implemented')
  }else{
    if(length(i) == 1){ ## they are just NAs
      if(only_slopes){
        .x = select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE) ## repeat the NAs
      }else{
        .x = python_like_select_name(i$par.fixed, "beta") ## repeat the NAs
      }
    }else{
      if(only_slopes){
        .x = select_slope_2(python_like_select_name(TMB::summary.sdreport(i)[,2], "beta"),v=F)
      }else{
        .x = python_like_select_name(TMB::summary.sdreport(i)[,2], "beta")
      }
    }
   .x
  }
}

createbarplot_object = function(fle, slotname='count_matrices_all', pre_path_funs=NULL){
  .x = readRDS(fle)
  if(is.null(pre_path_funs)) pre_path_funs <- "../../../CDA_in_Cancer/code/"
  lapply(.x, createbarplot_ROOSigs, slot=slotname, pre_path = pre_path_funs)
}


simulate_from_DM_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")[x_matrix[,2]+1])
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

simulate_from_M_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = softmax(cbind(logRmat, 0))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
  }
  return(sim_thetas)
}

simulate_from_correlated_binom = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, return_logratios=F){
  d = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, d)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    if(return_logratios){
      return(logRmat)
    }else{
      ## return first probability
      return(apply(logRmat, 2, function(i) exp(i)/(1+exp(i))))
    }
  }else{
    stop('Not implemented yet')
    # sim_thetas = sapply(1:d,
    #                     function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
    #                              give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    return(sim_thetas)
  }
}


give_confidence_interval = function(vec_est, vec_stderr){
  sapply(1:length(vec_est), function(i) c(vec_est[i]-1.96*vec_stderr[i],vec_est[i]+1.96*vec_stderr[i]) )
}

give_params_in_CI = function(vec_est, vec_stderr, vec_true){
  ci = give_confidence_interval(vec_est, vec_stderr)
  return(sapply(1:length(vec_est), function(i) (ci[1,i] < vec_true[i]) & (ci[2,i] > vec_true[i]) ))
}

wrapper_run_ttest_ilr = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 = compositions::ilr(props[[1]]), x2 = compositions::ilr(props[[2]]))$pvalue)
}
wrapper_run_ttest_props = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 =props[[1]][,-1], x2 = props[[2]][,-1])$pvalue)
}


summarise_DA_detection = function(true, predicted){
  warning('11 August 2021: there was a problem in how FP, TP, etc. were calculated')
  require(ROCR)
  ## remove NAs
  which_na = which(is.na(predicted))
  if(length(which_na) > 0){ ## some NA
    true = true[-which_na]
    predicted = predicted[-which_na]
  }
  
  FPs = sum(!true & predicted)
  FPR = FPs/sum(!true)
  TPs = sum(true & predicted)
  TPR = TPs/sum(true)
  TNs = sum(!true & !predicted)
  TNR = TNs/sum(!true)
  FNs = sum(true & !predicted)
  FNR = FNs/sum(!true)
  Accuracy=(TPs + TNs)/length(true)
  if(sum(true) == 0 | sum(!true) == 0){
    WeightedAccuracy = NA
  }else{
    WeightedAccuracy=mean(c(TPs/sum(true), TNs/sum(!true)))
  }
  total_pos = sum(true | predicted)
  Power = TPs/total_pos
  Sensitivity = TPs / (TPs + FNs)
  Specificity = TNs / (TNs + FPs)
  Recall=Sensitivity
  Precision=TPs/(TPs + FPs)
  pred <- (try(ROCR::prediction(as.numeric(true), as.numeric(predicted))))
  if(typeof(pred) == 'S4'){
    AUC = as.numeric(try(ROCR::performance(pred, "auc")@y.values[[1]]))
  }else{
    AUC = NA
  }
  return(c(FPR=FPR, TPR=TPR, Power=Power, AUC=AUC, Specificity=Specificity,
           Sensitivity=Sensitivity, Recall=Recall, Precision=Precision,
           Accuracy=Accuracy, WeightedAccuracy=WeightedAccuracy))
}

fill_covariance_matrix = function(arg_d, arg_entries_var, arg_entries_cov, verbose=T){
  if(verbose) warning('This function had been incorrect until now (30 july 2021)')
  .sigma <- give_UNSTRUCTURED_CORR_t_matrix(vec = arg_entries_cov, dim_mat = arg_d)
  # .sigma = matrix(NA, arg_d, arg_d)
  diag(.sigma) = arg_entries_var
  # .sigma[unlist(sapply(1:(arg_d-1), function(i) (i-1)*arg_d + (i+1):arg_d ))] = arg_entries_cov
  # .sigma[unlist(sapply(1:(arg_d-1), function(i) (i) + ((i):(arg_d-1))*arg_d))] = arg_entries_cov
  return(.sigma)
}

give_subset_sigs = function(sig_obj, sigs_to_remove){
  
  if(typedata %in% c("nucleotidesubstitution1", "nucleotidesubstitution3", "simulation")){
    slot_name = "count_matrices_all"
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(sig_obj,"count_matrices_active")[[1]]) | (length(attr(sig_obj,"count_matrices_active")[[1]]) == 0)){
      ## no active signatures
      slot_name = "count_matrices_all"
    }else{
      slot_name = "count_matrices_active"
    }
  }
  sig_obj_slot = attr(sig_obj,slot_name)
  sig_obj_slot = lapply(sig_obj_slot, function(i) i[, !(colnames(i) %in% sigs_to_remove)])
  # sig_obj$Y = sig_obj$Y[, !(colnames(sig_obj$Y) %in% sigs_to_remove)]
  attr(sig_obj,slot_name) = sig_obj_slot
  return(sig_obj)
}

give_amalgamated_exposures = function(sig_obj, list_groupings){
  
  if(typedata %in% c("nucleotidesubstitution1", "nucleotidesubstitution3", "simulation")){
    slot_name = "count_matrices_all"
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(sig_obj,"count_matrices_active")[[1]]) | (length(attr(sig_obj,"count_matrices_active")[[1]]) == 0)){
      ## no active signatures
      slot_name = "count_matrices_all"
    }else{
      slot_name = "count_matrices_active"
    }
  }
  sig_obj_slot = attr(sig_obj,slot_name)
  sig_obj_slot = lapply(sig_obj_slot, function(i){
    new_mat = sapply(list_groupings, function(j){
      grouped_exp <- i[,colnames(i) %in% j]
      if(!is.null(ncol(grouped_exp))){
        rowSums(grouped_exp)
      }else{
        grouped_exp
      }
      })
  colnames(new_mat) = sapply(list_groupings, function(i){if(length(i)>1){paste0(i[1], '+')}else{i}})
  return(new_mat)
  })
  attr(sig_obj,slot_name) = sig_obj_slot
  return(sig_obj)
}

give_amalgamated_exposures_TMBobj = function(sig_obj, list_groupings){
  sig_obj$Y = sapply(list_groupings, function(j){
    grouped_exp <- sig_obj$Y[,colnames(sig_obj$Y) %in% j]
    if(!is.null(ncol(grouped_exp))){
      rowSums(grouped_exp)
    }else{
      grouped_exp
    }
  })
  colnames(sig_obj$Y) = sapply(list_groupings, function(i){if(length(i)>1){paste0(i[1], '+')}else{i}})
  return(sig_obj)
}
  

give_subset_sigs_TMBobj = function(sig_obj, sigs_to_remove){
  sig_obj$Y = sig_obj$Y[,!(colnames(sig_obj$Y) %in% sigs_to_remove)]
  keep_obs = rowSums(sig_obj$Y) > 0
  sig_obj$Y = sig_obj$Y[keep_obs,]
  sig_obj$x = sig_obj$x[keep_obs,]
  sig_obj$z = sig_obj$z[keep_obs,]
  return(sig_obj)
}

normalise_rw <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise row-wise
    sweep(x, 1, rowSums(x), '/')
  }
}

normalise_cl <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise col-wise
    t(sweep(x, 2, colSums(x), '/'))
  }
}

give_barplot = function(ct, typedata, simulation=F, title='', legend_on=F, ...){
  require(gridExtra)
  obj = load_PCAWG(ct, typedata, simulation)
  a <- createBarplot(obj$Y[obj$x[,2] == 0,], remove_labels = T, ...)+ggtitle('Early raw')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  b <- createBarplot(obj$Y[obj$x[,2] == 1,], remove_labels = T, ...)+ggtitle('Late raw')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  c <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 0,]), remove_labels = T, ...)+ggtitle('Early normalised')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  d <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 1,]), remove_labels = T, ...)+ggtitle('Late normalised')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  
  if(!legend_on){
    a <- a+guides(fill=F)
    b <- b+guides(fill=F)
    c <- c+guides(fill=F)
    d <- d+guides(fill=F)
  }
  grid.arrange(a, b, c, d, top=title)
}

# wrapper_run_TMB_debug = function(object, model = "fullRE_DM", return_report=F, iter.max=150, init_log_lambda = 2, idx_cov_to_fill){
#   dim(object$Y)
#   # sort_columns=T
#   smart_init_vals=T
#   
#   data = object
#   
#   # if(sort_columns){
#   #   data$Y = data$Y[,order(colSums(data$Y), decreasing = F)]
#   # }
#   
#   data$Y = matrix(data$Y, nrow=nrow(data$Y))
#   data$x = (matrix(data$x, ncol=2))
#   
#   d <- ncol(data$Y)
#   n <- ncol(data$z) ## number of INDIVIDUALS, not samples
#   
#   data$num_individuals = n
#   data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#   
#   if(smart_init_vals){
#     require(nnet)
#     .x_multinom = multinom(data$Y ~ data$x[,2])
#     beta_init = t(coef(.x_multinom))
#   }else{
#     beta_init = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                         nrow = 2, byrow=TRUE))
#   }
#   
#   parameters <- list(
#     beta = beta_init,
#     u_large = matrix(runif(min = -0.3, max = 0.3, n = (d-1)*n), nrow=n),
#     logs_sd_RE=rep(1, d-1),
#     # cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#     cov_par_RE = runif(min = -1, max = 1, n = ((d-1)*(d-1)-(d-1))/2),
#     log_lambda = matrix(c(init_log_lambda,init_log_lambda))
#   )
#   if(model == "fullRE_DM"){
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "diagRE_DM"){
#     parameters$cov_par_RE = NULL
#     obj <- MakeADFun(data, parameters, DLL="diagRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "fullREDMsinglelambda"){
#     parameters$log_lambda = 2
#     obj <- MakeADFun(data, parameters, DLL="fullRE_dirichletmultinomial_single_lambda", random = "u_large")
#   }else if(model == "diagREDMsinglelambda"){
#     parameters$cov_par_RE = NULL
#     parameters$log_lambda = 2
#     obj <- MakeADFun(data, parameters, DLL="diagRE_dirichletmultinomial_single_lambda", random = "u_large")
#   }else if(model == "sparseRE_DM"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomial", random = "u_large")
#     }else{
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomial_single", random = "u_large")
#     }
#   }else if(model == "sparseRE_DMSL"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda", random = "u_large")
#     }else{
#       stop()
#     }
#   }else if(model == "sparseRE_DMSL2"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda2", random = "u_large")
#     }else{
#       stop()
#     }
#   }else if(model == "sparseRE_DMSL2"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda2", random = "u_large")
#     }else{
#       stop()
#     }
#   } else{
#     stop("Incorrect <model>")
#   }
# 
#   opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max, trace=T)
#   
#   if(return_report)  return(sdreport(obj))
#   return(opt)
# }

is_slope = function(v){
  bool_isbetaslope = rep(F, length(v))
  bool_isbetaslope[v == "beta"] = T
  bool_isbetaslope[which(bool_isbetaslope)] = c(F,T)
  bool_isbetaslope
}

give_barplot_from_obj <- function(obj, legend_on=F, legend_bottom=F, nrow_plot=2, only_normalised=F, title=NULL, ...){
  a <- createBarplot(obj$Y[obj$x[,2] == 0,], remove_labels = T, ...)+ggtitle('Early raw')
  b <- createBarplot(obj$Y[obj$x[,2] == 1,], remove_labels = T, ...)+ggtitle('Late raw')
  c <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 0,]), remove_labels = T, ...)+ggtitle('Early normalised')
  d <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 1,]), remove_labels = T, ...)+ggtitle('Late normalised')
  if(!legend_on){
    a <- a+guides(fill=F)
    b <- b+guides(fill=F)
    c <- c+guides(fill=F)
    d <- d+guides(fill=F)
  }
  if(legend_bottom){
    a <- a+theme(legend.position='bottom', legend.text = element_text(size = 6))
    b <- b+theme(legend.position='bottom', legend.text = element_text(size = 6))
    c <- c+theme(legend.position='bottom', legend.text = element_text(size = 6))
    d <- d+theme(legend.position='bottom', legend.text = element_text(size = 6))
  }
  if(only_normalised){
    grid.arrange(c, d, top=title, nrow=nrow_plot)
  }else{
    grid.arrange(a, b, c, d, top=title, nrow=nrow_plot)
  }
}

give_ranked_plot_simulation = function(tmb_fit_object, data_object, print_plot = T, nreps = 1, model, integer_overdispersion_param){
  
  if(model == 'M'){
    ## theta is always going to be the same. Only replicate the draws from the multinomial
    
    sim_theta = simulate_from_M_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                    x_matrix = data_object$x, z_matrix = data_object$z)
    
    sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
    
    if(nreps>1){
      sim_counts = replicate(nreps, t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) ), simplify = F)
    }
  }else if(model %in% c('DM', 'DMSL')){
    give_sim_data = function(){
      if(model == 'DM'){
        sim_theta = simulate_from_DM_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                         x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
      }else if(model == 'DMSL'){
        sim_theta = simulate_from_DMSL_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                         x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
      }
      sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
      return(sim_counts)
    }    
    if(nreps>1){
      sim_counts = replicate(nreps, give_sim_data(), simplify = F)
    }else{
      sim_counts = give_sim_data()
    }
  }else{
    stop('Specify a correct model')
  }
  
  stopifnot(all(dim(data_object$Y) == dim(sim_counts)))
  if(print_plot)  plot(sort(data_object$Y), sort(sim_counts))
  
  return(sim_counts)
  
}


simulate_from_M_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = softmax(cbind(logRmat, 0))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
  }
  return(sim_thetas)
}

simulate_from_DM_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")[x_matrix[,2]+1])
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

simulate_from_DMSL_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = rep(integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")), nrow(x_matrix))
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    ## univariate RE
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

give_interval_plots = function(df_rank, data_object, loglog=F){
  xx = melt(df_rank, id.vars=c('sorted_value', 'rank_number'))
  xx_summary = xx %>% group_by(rank_number) %>% mutate(min_interval=quantile(sorted_value, probs = c(0.025)),
                                                       max_interval=quantile(sorted_value, probs = c(0.975)),
                                                       mean=mean(sorted_value))
  xx_summary = xx_summary[!(duplicated(xx_summary[,c('rank_number', 'min_interval', 'max_interval')])),c('rank_number', 'min_interval', 'max_interval', 'mean')]
  a = ggplot(cbind.data.frame(xx_summary, sorted_true=sort(data_object$Y)), aes(x=sorted_true, ymin= min_interval, ymax=max_interval))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+
    geom_ribbon(fill = "red", alpha=0.2)+
    geom_point(aes(x=sorted_true, y=mean), size=0.4)+
    geom_line(aes(x=sorted_true, y=mean))+
    labs(x='Observed ranked value', y='Mean simulated ranked value')
  if(loglog){
    a = a+ scale_y_continuous(trans = "log")+scale_x_continuous(trans = "log")
  }
  
  return(a)
}

give_interval_plots_2 = function(df_rank, data_object,loglog=F, title, theme_bw=F){
  xx = melt(df_rank, id.vars=c('sorted_value', 'rank_number'))
  xx_summary = xx %>% group_by(rank_number) %>% mutate(min_interval=quantile(sorted_value, probs = c(0.025)),
                                                       max_interval=quantile(sorted_value, probs = c(0.975)),
                                                       mean=mean(sorted_value))
  xx_summary = xx_summary[!(duplicated(xx_summary[,c('rank_number', 'min_interval', 'max_interval')])),c('rank_number', 'min_interval', 'max_interval', 'mean')]
  xx_summary <- cbind.data.frame(xx_summary, sorted_true=as.vector(data_object$Y))
  xx_summary[,'rank_true'] = order(xx_summary$sorted_true)
  xx_summary$col = apply(xx_summary, 1, function(i) (i['sorted_true'] > i['min_interval']) & (i['sorted_true'] < i['max_interval']) )
  a  = ggplot(xx_summary, aes(x=sorted_true,
                              ymin= min_interval, ymax=max_interval))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+
    geom_ribbon(fill = "red", alpha=0.2)+
    geom_point(aes(x=sorted_true, y=mean, col=col), size=0.4)+
    geom_line(aes(x=sorted_true, y=mean))+
    labs(x='Observed value', y='Mean simulated value')+
    ggtitle(title, subtitle=paste0(paste0(names(table(xx_summary$col)), ':', table(xx_summary$col)), collapse ='; '))+
    theme(legend.position = "bottom")
  
  if(theme_bw){
    a = a+theme_bw()+theme(legend.position = "bottom")
  }
  
  if(loglog){a=a+ scale_y_continuous(trans = "log")+scale_x_continuous(trans = "log") }
  
  return(a)
}

give_betas <- function(TMB_obj){
  matrix(python_like_select_name(TMB_obj$par.fixed, 'beta'), nrow=2)
}



plot_betas <- function(TMB_obj, names_cats=NULL, rotate_axis=T, theme_bw=T, remove_SBS=T, only_slope=F, return_df=F, plot=T,
                       line_zero=T, add_confint=F, return_plot=T, return_ggplot=F, title=NULL){
  if(typeof(TMB_obj) == 'character'){
    .summary_betas <- NA
    if(theme_bw){
      plt <- ggplot()+theme_bw()
      if(plot) print(plt)
    }else{
      plt <- ggplot()
      if(plot) print(plt)
    }
  }else{
    .summary_betas <- summary(TMB_obj)
    .summary_betas <- cbind.data.frame(python_like_select_rownames(.summary_betas, 'beta'),
                                       type_beta=rep(c('Intercept', 'Slope')),
                                       LogR=rep(1:(nrow(python_like_select_rownames(.summary_betas, 'beta'))/2), each=2))
    if(only_slope){
      .summary_betas <- .summary_betas[.summary_betas$type_beta == 'Slope',]
    }
    
    if(!is.null(names_cats)){
      if(remove_SBS){
        names_cats <- gsub("SBS", "", names_cats) 
      }
      if(length(unique(.summary_betas$LogR)) != length(names_cats)){
        stop('Number of beta slope/intercept pairs should be the same as the length of the name of the categories')
      }
      .summary_betas$LogR = names_cats[.summary_betas$LogR]
    }
    plt <- ggplot(.summary_betas, aes(x=LogR, y=`Estimate`))
    
    if(line_zero) plt <- plt + geom_hline(yintercept = 0, lty='dashed', col='blue')
      
    plt <- plt +
      geom_point()+
      geom_errorbar(aes(ymin=`Estimate`-`Std. Error`, ymax=`Estimate`+`Std. Error`), width=.1)+
      ggtitle('Slopes')+facet_wrap(.~type_beta, scales = "free")
    
    if(theme_bw){
      plt <- plt + theme_bw()
    }
    
    if(add_confint){
      confints <- cbind(.summary_betas, confint=t(give_confidence_interval(.summary_betas[,'Estimate'], .summary_betas[,'Std. Error'])))
      plt <- plt+
        geom_errorbar(data = confints, aes(ymin=confint.1, ymax=confint.2), width=.1 ,col='blue', alpha=0.6)

    }
    
    if(rotate_axis){
      plt <- plt + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    
    if(!is.null(title)){
      plt <- plt + ggtitle(title)
    }
    
    if(!TMB_obj$pdHess){
      plt <- plt + annotate("text", x = 0, y=.5, label="not PD")+geom_point(col='red')
      if(plot) print(plt)
    }else{
      if(plot) print(plt)
    }
  }
  
  if(return_df){
    .summary_betas
  }else{
    if(return_plot & return_df){stop('<return_plot=T> and <return_df=T> are incompatible')}
    plot_list <- list(plt)
    class(plot_list) <- c("quiet_list", class(plot_list))
    if(return_plot){
      return(cowplot::as_grob(plt))
    }else if(return_ggplot){
      return(plt)
    }
  }
}


plot_lambdas <- function(TMB_obj, return_df=F, plot=T){
  
  lambdas_df <- python_like_select_rownames(summary(TMB_obj), 'log_lambda')
  if(!is.null(dim(lambdas_df))){
    .summary_lambda <- cbind.data.frame(data.frame(lambdas_df),
                                        name=c('Lambda 1', 'Lambda 2'))
  }else{
    .summary_lambda <- cbind.data.frame(t(lambdas_df), name='Lambda 1')
    colnames(.summary_lambda) <- make.names(colnames(.summary_lambda))
  }
  
  plt <- (ggplot(.summary_lambda, aes(x=name, y=`Estimate`))+
            geom_point()+
            geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw())
  if(plot){
    print(plt)
  }
  
  if(return_df){
    return(.summary_lambda)
  }else{
    return(plt)
  }
}

plot_estimates_TMB <- function(TMB_obj, parameter_name, return_df=F, plot=T, verbatim=T){
  if(verbatim) cat('Consider the other functions <plot_betas> and <plot_lambdas>\n')
  
  
  parameter_df <- python_like_select_rownames(summary(TMB_obj), parameter_name)
  if(length(parameter_df) == 0){
    ## there is no parameter
    return(cbind.data.frame(Estimate=NA, `Std..Error`=NA, name=NA))
  }else{
    if(!is.null(dim(parameter_df))){
      .summary_param <- cbind.data.frame(data.frame(parameter_df),
                                          name=paste0(parameter_name, 1:nrow(parameter_df)))
    }else{
      ## single parameter
      .summary_param <- cbind.data.frame(t(parameter_df), name=paste0(parameter_name, '1'))
      colnames(.summary_param) <- make.names(colnames(.summary_param))
    }
  
    plt <- (ggplot(.summary_param, aes(x=name, y=`Estimate`))+
              geom_point()+
              geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw())
    if(plot){
      print(plt)
    }
  
    if(return_df){
      return(.summary_param)
    }else{
      return(plt)
    }
  }
}

give_length_cov <- function(dim1_covmat){
  ((dim1_covmat**2) - dim1_covmat )/2
}


give_sim_from_estimates <- function(ct, typedata = "signatures", sigs_to_remove="", model="sparseRE_DM",
                                    bool_nonexo=TRUE, bool_give_PCA, sig_of_interest='SBS8',
                                    path_to_data= "../../../data/", tmb_object=NULL, obj_data=NULL,
                                    nrow_pca_plot=2){
  
  if(is.null(tmb_object)){
    if(model == "fullRE_M"){
      if(!bool_nonexo)    list_estimates <- fullRE_M
      if(bool_nonexo)    list_estimates <- fullRE_M_nonexo
    }else if(model == "fullRE_DM"){
      if(!bool_nonexo)    list_estimates <- fullRE_DMSL
      if(bool_nonexo)    list_estimates <- fullRE_DMSL_nonexo
    }else if(model == "sparseRE_DM"){
      if(bool_nonexo)    list_estimates <- sparseRE_DMSL_nonexo
    }
  }else{
    list_estimates <- list()
    list_estimates[[ct]] <- tmb_object
  }
  
  if(is.null(obj_data)){
    warning('WARNING! Here I am sorting the columns of TMB, I should not in all cases. Specify <obj_data> if needed')
    obj_data <- sort_columns_TMB(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = typedata, path_to_data =path_to_data),
                                                       sigs_to_remove = sigs_to_remove))
  }
  dmin1 <- ncol(obj_data$Y)-1
  cov_vec = rep(0, (dmin1**2-dmin1)/2)
  
  if(model %in% c("fullRE_M", "fullRE_DM", "fullRE_DMSL")){
    cov_vec = python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_par_RE')
    ###**** I AM NOT SURE ABOUT THIS BIT BELOW! ARE THEY SD OR VAR???*****###
    ### implementing them as though they were sd ###
    var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))**2
    var_vec_v2 = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))
  }else if(model == "sparseRE_DM"){
    cov_vec[as.numeric(strsplit(subset_sigs_sparse_cov_idx_nonexo[subset_sigs_sparse_cov_idx_nonexo$V1 == ct,"V2"], ',')[[1]])] = python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_RE_part')
    var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))**2
    var_vec_v2 = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))
  }else if(model %in% c("diagRE_M", "diagRE_DM", "diagRE_DMSL")){
    cov_vec = rep(0, give_length_cov(length(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))))
    var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))**2
    var_vec_v2 = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))
  }else{
    stop('Specify correct <model>')
  }
  
  cov_mat <- fill_covariance_matrix(arg_d = dmin1,
                                    arg_entries_var = var_vec,
                                    arg_entries_cov = cov_vec)
  cov_mat_v2 <- fill_covariance_matrix(arg_d = dmin1,
                                       arg_entries_var = var_vec_v2,
                                       arg_entries_cov = cov_vec)
  # cov_matb <- fill_covariance_matrix(arg_d = dmin1,
  #                                   arg_entries_var = var_vec**2,
  #                                   arg_entries_cov = cov_vec)
  
  beta_mat = matrix(python_like_select_name(list_estimates[[ct]]$par.fixed, 'beta'), nrow=2)
  
  n_sim = 1000
  x_sim = cbind(1, rep(c(0,1), n_sim))
  u_sim = mvtnorm::rmvnorm(n = n_sim, mean = rep(0,dmin1), sigma = cov_mat)
  
  theta = x_sim %*% beta_mat + (give_z_matrix(n_sim*2)) %*% u_sim
  
  if(model %in% c('sparseRE_DM', 'fullRE_DM', 'fullRE_DMSL', 'diagRE_DMSL')){
    alpha = softmax(cbind(theta, 0))*exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'log_lambda'))
  }else if(model %in% c("fullRE_M")){
    alpha = softmax(cbind(theta, 0))
  }else{
    stop('Check softmax step')
  }
  
  if(model %in% c('sparseRE_DM', 'fullRE_DM', 'fullRE_DMSL', 'diagRE_DMSL')){
    probs = t(apply(alpha, 1, MCMCpack::rdirichlet, n=1))
  }else if(model %in% c("fullRE_M")){
    probs = alpha
  }else{
      stop('Check probabilities step')
    }
  
  probs_obs = normalise_rw(obj_data$Y)
  all_probs = rbind(probs_obs, probs)
  
  if(bool_give_PCA){
    pca <- prcomp(all_probs)
    
    if(! (sig_of_interest %in% colnames(all_probs))){stop('Specify a <sig_of_interest> present in the dataset')}
    
    df_pca <- cbind.data.frame(pca=pca$x[,1:2], col=c(rep('Observed', nrow(probs_obs)), rep('Simulated',nrow(probs))),
                               sig_of_interest=all_probs[,sig_of_interest],
                               group=c(rep(c('early','late')[obj_data$x[,2]+1]),
                                       rep(c('early','late'), n_sim)))
    return(list(df_pca, ggplot(df_pca, aes(x=pca.PC1, y=pca.PC2, col=sig_of_interest))+labs(col=sig_of_interest)+
                  geom_point(alpha=0.7)+facet_wrap(.~interaction(col,group),nrow=nrow_pca_plot)))
  }else{
    return(all_probs)
  }
  
}

vector_cats_to_logR <- function(i){paste0(i[-length(i)], '/', i[length(i)])}

non_duplicated_rows <- function(i){
  rownames(i)[duplicated(rownames(i))] = paste0(rownames(i)[duplicated(rownames(i))], "_2")
  i
}

give_all_correlations <- function(j){
  outer(1:nrow(j), 1:nrow(j), Vectorize(function(i,k){cor(j[i,], j[k,])}))
}

give_perturbation_TMBobj = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_perturbation_ROOSigs_alt

  if(addone){
    exposures_cancertype_obj$Y = exposures_cancertype_obj$Y + 1
  }
  
  ## for each individual, compute the perturbation
  apply(exposures_cancertype_obj$z, 2, function(i){
    ## if there are two individuals
    if(sum(i) == 2){
      exposures_cancertype_obj$Y[which(i == 1)[1],] / exposures_cancertype_obj$Y[which(i == 1)[2],]
    }
  })
  
}

give_totalperturbation_TMBobj = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_total_perturbation
  
  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  mean(apply(pert, 2, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2))))
  
}

give_totalperturbation_TMBobj_sigaverage = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_total_perturbation
  
  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  mean(apply(pert, 1, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2))))
  
}

give_weightedtotalperturbation_TMBobj = function(exposures_cancertype_obj, addone=F){

  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  perts <- apply(pert, 2, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2)))
  weights <- as.vector(apply(exposures_cancertype_obj$z, 2, function(i){
    ## if there are two individuals
    if(sum(i) == 2){
    log(sum(rowSums(exposures_cancertype_obj$Y[which(i == 1),])))
    }
  }))
  if(length(perts) != length(weights)){
    stop('Weights and perturbations do not have the same length')
  }
  weights <- weights/sum(weights) # normalise
  sum(perts * weights)/nrow(exposures_cancertype_obj$Y)
  
}

wrapper_run_HMP_Xdc.sevsample <- function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  return(HMP::Xdc.sevsample(x)$`p value`)
}


wrapper_run_HMP_Xmcupo.sevsample <- function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  return(HMP::Xmcupo.sevsample(x)$`p value`)
}

plot_ternary <- function(x, legend_on=T, plot_points=T, ...){
  require(Ternary)
  
  if(ncol(x) != 3){stop('Number of columns must be three. Create a subcomposition or amalgamation if needed.')}
  TernaryPlot(atip = colnames(x)[1], btip = colnames(x)[2], ctip = colnames(x)[3],
              grid.lines = 0, grid.col = NULL, ...)
  dens <- TernaryDensity(x, resolution = 10L)
  
  cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                     seq(from = 0, to = 47, by=1))
  if(legend_on){
    legend(x=-0.4,y=1.08,
           fill = cls_legend[1,][c(T,F,F,F,F)],
           legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
           y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
  }
  ColourTernary(dens)
  if(plot_points)  TernaryPoints(x, col = 'red', pch = '.', cex=5)
  TernaryDensityContour(x, resolution = 30L)
}

split_matrix_in_half <- function(x){
  list(x[1:(nrow(x)/2),],
       x[(1+(nrow(x)/2)):nrow(x),])
}

