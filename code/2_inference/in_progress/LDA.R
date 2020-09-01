## Infence using stan_logistic_normal_multinomial_ME
## with both random effects and fixed effects

rm(list = ls())
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Sys.setenv(LANG='en')

# set.seed(123)
# set.seed(223)

library(rstan)
#library(uuid)
#library(optparse)
#library(dplyr)
#library(TMB)
#library(rgl) ## for 3D
#library(gridExtra)
source("~/Documents/PhD/GlobalDA/code/2_inference/helper/helper_DA_stan.R")
source("~/Documents/PhD/GlobalDA/code/1_create_ROO/roo_functions.R")
source("~/Documents/PhD/GlobalDA/code/1_create_ROO/helper_1_create_ROO.R")

simulation=FALSE
if(simulation){
  opt=list()
  opt$input = NULL
  opt$d = 4
  opt$n = 10
  opt$nlambda = 50
  opt$beta_gamma_shape = 7
  opt$infile='~/Documents/PhD/GlobalDA/data/assessing_models_simulation/datasets/20200625_10_100_3_0_dataset.RDS'
  opt$output = NULL
  opt$niterations = 10000
  opt$model = 'LDA_FE'
  Nits = opt$niterations
  
  # objects_sigs_per_CT_features <- readRDS(opt$infile)
  
  
  ## (A) simulate from the model with different parameters
  d = opt$d # opt$Nk ## number of signatures
  n = opt$n #opt$Ns ## number of samples
  beta_gamma_shape = opt$beta_gamma_shape  ##opt$hyperparam_shape ## shape parameter for the beta
  sd_RE = 0.3 ## standard deviation for random effects
  lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
  Nm_lambda = opt$nlambda ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)
  
  ## Group effects
  ## covariate matrix
  X_sim = matrix(NA, nrow=2, ncol=2*n)
  ## the samples are split into two groups
  X_sim[1,] = 1
  X_sim[2,] = rep(c(0,1), each=n)
  beta = matrix(0, nrow=2, ncol=d-1)
  beta[2,] = rgamma(n = d-1, shape = beta_gamma_shape, rate = beta_gamma_shape) ## for the coefficients
  
  ## Random effects
  Z_sim0 = matrix(0, nrow=n, ncol=n)
  diag(Z_sim0) = 1
  Z_sim = t(rbind(Z_sim0, Z_sim0))
  u = matrix(NA, nrow=n, ncol=1)
  u[,1] = rnorm(n = n, mean = 0, sd = sd_RE)
  
  ## lambdas: overdispersion
  lambdas = c(rep(lambda[1], n), rep(lambda[2], n))
  
  ## create alpha
  alphabar = softmax_mat( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )
  
  emit_DM = FALSE
  if(emit_DM){
    alpha = alphabar * lambdas
    
    # image(t(alpha), main='Alphas')
    
    ## Create the counts
    Nm = rpois(n*2, lambda = Nm_lambda) # number of mutations per sample \in N^{2*n}
    W = matrix(NA, nrow = 2*n, ncol = d)
    for(l in 1:(2*n)){
      W[l,] = HMP::Dirichlet.multinomial(Nrs = Nm[l], shape = alpha[l,])
    }
  }else{
    theta = alphabar
    W = t(apply(theta, 1, rmultinom, n = 1, size = Nm_lambda))
    Nm= Nm_lambda
  }
  
#  image(t(W))
  
  ## Save as object so that we can perform the inference
  objects_counts <- new("exposures_cancertype",
                        cancer_type="simulated data",
                        type_classification = "simulated two group",
                        number_categories = 2,
                        id_categories = c('sim1', 'sim2'),
                        active_signatures = "absent; simulation",
                        count_matrices_all = list(give_dummy_row_names(W[1:n,]), give_dummy_row_names(W[(n+1):(2*n),])),
                        count_matrices_active = list(list(), list()),
                        sample_names = rownames(give_dummy_row_names(W[1:n,])),
                        modification = "none",
                        is_null_active = TRUE,
                        is_empty="Non-empty"
  )
  
  objects_sigs_per_CT_features = list(objects_counts=objects_counts, d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=Nm_lambda,
                                      X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = Nm, W = W)
  
  d = objects_sigs_per_CT_features$d ## number of signatures or features
  n = objects_sigs_per_CT_features$n ## number of samples
  X = objects_sigs_per_CT_features$X
  Z = objects_sigs_per_CT_features$Z
  objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features$objects_counts,"count_matrices_all")
  W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]])
}else{
    ## read PCAWG data
    opt = list()
    opt$infile = "/Users/morril01/Documents/PhD/GlobalDA/data/roo/Biliary-AdenoCA_signatures_ROO.RDS"
    opt$cancer_type = strsplit(basename(opt$infile), '_')[[1]][1]
    opt$typedata = strsplit(basename(opt$infile), '_')[[1]][2]
    opt$output = NULL
    opt$niterations = 200
    opt$model = "LDA"
    
    Nits = opt$niterations
    
    fle <- opt$infile
    if(file.exists(fle)){
      objects_sigs_per_CT_features <- readRDS(fle)
    }else{
      cat('Object for cancer type', opt$cancertype, 'was not found')
    }
    
    if(opt$typedata == "nucleotidesubstitution1"){
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else if(opt$typedata == "nucleotidesubstitution3"){
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else if(opt$typedata == "signatures"){
      if(is.null(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]])){
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
    
  }

model_file_name = paste0(opt$model, ".stan")
rstan_options(auto_write = TRUE)
stanc(model_file_name)


flatten_matrix = function(x){
  ## give a vector with ids from a matrix. it is performed column-wise
  as.vector(unlist(sapply(1:ncol(x), function(i) unlist(sapply(x[,i], function(j) rep(i, j))))))
}

give_sample_idx_matrix = function(x){
  unlist(sapply(1:ncol(x), function(i) rep(1:nrow(x), x[,i])))
}

stan_data = list(M=nrow(W),
                 w = flatten_matrix(W),
                 doc = give_sample_idx_matrix(W),
                 V = ncol(W),
                 K=2,
                 N = sum(W),
                 alpha = rep(1, 2),
                 beta = rep(1, ncol(W)),
                 # p=2,
                 # w = W,
                 x = X)
                 # Z = Z)

# params = c('beta', 'u')

fit_stan <- stan(file = model_file_name, data = stan_data,
                 iter = Nits, chains = 4, cores = 4, thin = 1, control = list(stepsize=1, adapt_delta=0.98, max_treedepth=12))
# Warning messages:
#   1: There were 12062 divergent transitions after warmup. Increasing adapt_delta above 0.98 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 7424 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 12. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems

#save.image(file = "~/Desktop/image_LDA_FE_PCAWG.RData")
save.image(file = "/home/lm687/simple_LDA_image.RData") 
# save.image(file = "~/Desktop/image_LDA_FE.RData")

# fit_stan2 = (extract(fit_stan))
# 
# ## all pairs
# pairs(fit_stan, pars = names(fit_stan))
# ## all beta coefficients
# pairs(fit_stan, pars = names(fit_stan)[grep('beta', names(fit_stan))])
# ## beta coefficients for intersect
# pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[1,', names(fit_stan))], text.panel = "Coefficients for the intercept")
# ## beta coefficients for slope
# pairs(fit_stan, pars = names(fit_stan)[grep('beta\\[2,', names(fit_stan))], text.panel = "Coefficients for the slope")
# 
# dim(fit_stan2$beta)
# mean(fit_stan2$beta[,1,1])
# mean(fit_stan2$beta[,1,2])
# mean(fit_stan2$beta[,1,3])
# 
# pairs(cbind(fit_stan2$beta[,1,1],fit_stan2$beta[,1,2], fit_stan2$beta[,1,3]))
# 
# test <- invisible(nnet::multinom(W ~as.factor(X[2,])))
# 
# beta
# t(sapply(1:2, function(i) colMeans(fit_stan2$beta[,i,] )))
# t(coefficients(test))
# 
# pairs(cbind(true=beta %>% as.vector,
#             stan_mean=t(sapply(1:2, function(i) colMeans(fit_stan2$beta[,i,] ))) %>% as.vector,
#             nnet=t(coefficients(test)) %>% as.vector))
# 
# ## 3d of LP
# open3d() %>% invisible()
# rgl::plot3d(fit_stan2$beta[,2,][,1], fit_stan2$beta[,2,][,2], fit_stan2$lp__, xlab = "w", ylab = "k", zlab = "lp__")
# summary(fit_stan)$summary
# 
# ## compare posterior and prior
# source("~/Documents/PhD/other_repos/b_tape/Vias_Brenton/copy_number_analysis_organoids/helper_functions.R")
# plts_posterior_prior = sapply(1:dim(fit_stan2$beta)[3], function(idx_3){
#   lapply(1:dim(fit_stan2$beta)[2], function(idx_2){
#     give_joint_histogram(some_list = list(posterior=fit_stan2$beta[,idx_2,1],
#                                           prior=runif(1000, -5, 5)))+ggtitle(paste0('beta ', idx_2, ' ', idx_3))
#   })
# })
# 
# plts_posterior_prior = c(plts_posterior_prior,
#                          lapply(1:dim(fit_stan2$u)[2], function(idx_2){
#                            (give_joint_histogram(some_list = list(posterior=fit_stan2$u[,idx_2],
#                                                                   prior=rnorm(1000, mean = 0, sd = rgamma(1, 5,5))))+ggtitle(paste0( 'u', idx_2)))
#                          }))
# 
# library(gridExtra)
# do.call(grid.arrange, plts_posterior_prior)
# 
# 
# compile("~/Desktop/TMB/mm_multinomial/FE_multinomial.cpp", "-std=gnu++17")
# dyn.load(dynlib("~/Desktop/TMB/mm_multinomial/FE_multinomial"))
# data <- list(Y = W, n=n, d=d, x=t(X))
# parameters <- list(
#   beta = array(rbind(rep(1,d-1), rep(1,d-1)), dim = c(2,d-1))
# )
# obj <- MakeADFun(data, parameters, DLL="FE_multinomial")
# sdreport(obj)
# 
# # betas_estimate = rbind(opt$par[grep('beta_intersect', names(opt$par))],
# #                        c(opt$par[grep('beta_slope', names(opt$par))]))
# betas_estimate = matrix(opt$par, nrow=2)
# betas_estimate
