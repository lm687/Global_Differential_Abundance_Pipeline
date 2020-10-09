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

load_PCAWG = function(ct, typedata, simulation=FALSE){
  if(simulation){
    fle = ct
  }else{
    fle = paste0("../../data/roo/", ct, '_', typedata, "_ROO.RDS" )
  }
  objects_sigs_per_CT_features <- readRDS(fle)
  
  if(simulation){
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
  }else if(typedata == "signatures"){
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
  }
  
  d = ncol(objects_sigs_per_CT_features[[1]]) ## number of signatures or features
  n = nrow(objects_sigs_per_CT_features[[1]]) ## number of samples
  
  print(n)
  
  # covariate matrix
  X = matrix(NA, nrow=2, ncol=2*n)
  X[1,] = 1
  X[2,] = rep(c(0,1), each=n)
  
  Z0 = matrix(0, nrow=n, ncol=n)
  diag(Z0) = 1
  Z = t(rbind(Z0, Z0))
  
  ## The counts
  W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])
  
  return(list(x=t(X), z=t(Z), Y=W))
}

wrapper_run_TMB = function(ct, typedata, model, simulation=FALSE, allow_new_LNM=FALSE){
  cat(ct)
  cat(typedata)
  data = load_PCAWG(ct, typedata, simulation)
  if(length(data) == 1){
    if(is.na(data)){
      return(warning('RDS object is NA'))
    }
  }
  data$Y = matrix(data$Y, nrow=nrow(data$Y))
  data$x = (matrix(data$x, ncol=2))
  
  d <- ncol(data$Y)
  n <- ncol(data$z) ## number of INDIVIDUALS, not samples
  
  if(model == "M"){
    data$num_individuals = n
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_random_effects = matrix(rep(1, n)),
      logSigma_RE=1
    )
    obj <- MakeADFun(data, parameters, DLL="ME_multinomial", random = "u_random_effects")
  }else if(model == "fullRE_M"){
    data$num_individuals = n
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial", random = "u_large")
  }else if(model == "diagRE_M"){
    data$num_individuals = n
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_multinomial", random = "u_large")
  }else if(model == "LNM"){
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_random_effects = matrix(rep(1, n)),
      logSigma_RE=0,
      cov_par = rep(1, ((d-1)*(d-1)-(d-1))/2)
    )
    obj <- MakeADFun(data, parameters, DLL="ME_LNM", random = "u_random_effects")
  }else if(model == "fullRE_LNM"){
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logSigma_RE=0,
      cov_par = rep(1, ((d-1)*(d-1)-(d-1))/2),
      logs_sd_RE=rep(1, d-1)
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_LNM", random = "u_large")
  }else if(model == "DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_random_effects = matrix(rep(1, n)),
      logSigma_RE=1,
      log_lambda = matrix(c(2,2))
    )
    obj <- MakeADFun(data, parameters, DLL="ME_dirichletmultinomial", random = "u_random_effects")
  }else if(model == "fullRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
      log_lambda = matrix(c(2,2))
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial", random = "u_large")
  }else if(model == "fullRE_DM_singlelambda"){
    data$num_individuals = n

    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
      log_lambda = 1
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_singlelambda_dirichletmultinomial", random = "u_large")
  }else if(model == "diagRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      log_lambda = matrix(c(2,2))
    )
    obj <- MakeADFun(data, parameters, DLL="diagRE_ME_dirichletmultinomial", random = "u_large")
  }else if(model == "fullRE_DM_altpar"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- list(
      beta = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                     nrow = 2, byrow=TRUE)),
      u_large = matrix(rep(1, (d-1)*n), nrow=n),
      logs_sd_RE=rep(1, d-1),
      cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
      log_lambda = matrix(c(2,2))
    )
    obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial_altpar", random = "u_large")
  }else{
    stop('Specify correct <model>\n')
  }
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  # obj$he()    ## <-- Analytical hessian
  return(sdreport(obj))
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
  matrix(vec_RE, ncol=dmin1)
}

give_summary_per_sample = function(TMB_object){
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

give_summary_of_runs = function(vector_TMB_objects, long_return){
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
  chisqrt_stat = t(v) %*% solve(sigma) %*% v
  pchisq(q = chisqrt_stat, df = length(v), lower.tail = FALSE)
}

wald_TMB_wrapper = function(i, verbatim=TRUE){
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = select_slope_2(which(names(i$par.fixed) == "beta"), verbatim=verbatim)
    wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
  }
}


select_slope = function(i){
  stop('Deprecated due to error! use <select_slope_2> instead')
  i[(length(i)/2+1):(length(i))]
}

select_slope_2 = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    if(verbatim) warning('As per 27 August it seems clear that this version, and not <select_slope>, is correct')
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

