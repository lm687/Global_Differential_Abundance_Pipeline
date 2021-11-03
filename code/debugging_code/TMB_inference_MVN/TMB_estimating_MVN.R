##----------------------------------------------------------------##
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(MASS)
library(mvtnorm)
library(TMB)
library(matrixcalc)
source("../../../GlobalDA/code/2_inference_TMB/helper_TMB.R")

# TMB::compile("tmb_MVN.cpp", "-std=gnu++17")
# dyn.load(dynlib("tmb_MVN"))
TMB::compile("tmb_MVN_nolog.cpp", "-std=gnu++17")
dyn.load(dynlib("tmb_MVN_nolog"))
##----------------------------------------------------------------##

##----------------------------------------------------------------##
## Functions
give_UNSTRUCTURED_CORR_t_matrix = function(vec, dim_mat){
  ##' given the off-diagonal, non-redundant, entries of a symmetric matrix,
  ##' populate the matrix. Diagonal entries are doing to be set to 1
  # #https://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html
  m = matrix(1, nrow = dim_mat, ncol = dim_mat)
  ## fill in the order that TMB's UNSTRUCTURED_CORR_t saves the covariances
  m[unlist(sapply(2:nrow(m), function(rw) seq(from = rw,length.out = (rw-1), by = nrow(m) )))] = vec
  m[unlist(sapply(2:nrow(m), function(cl) seq(from = (((cl-1)*nrow(m))+1),length.out = (cl-1), by = 1 )))] = vec
  return(m)
}

fill_covariance_matrix = function(arg_d, arg_entries_var, arg_entries_cov){
  ##' given the off-diagonal, non-redundant, entries of a symmetric matrix,
  ##' as well as the diagonal entries, populate the covariance matrix
  .sigma <- give_UNSTRUCTURED_CORR_t_matrix(vec = arg_entries_cov, dim_mat = arg_d)
  diag(.sigma) = arg_entries_var
  return(.sigma)
}

python_like_select_name = function(vector, grep_substring){
  ## subset vector with entries that match name
  vector[grepl(pattern = grep_substring, x = names(vector))]
}
##----------------------------------------------------------------##


##----------------------------------------------------------------##
## simulating data from a four-dimensional MVN
d <- 4
example_mat <- matrix(runif(d^2)*2-1, ncol=d)
## sigma_true is our covariance matrix
Sigma_true <- t(example_mat) %*% example_mat
matrixcalc::is.positive.definite(Sigma_true)
## the true mean is a vector of zeros
Y = mvtnorm::rmvnorm(n = 5000, mean = rep(0,d), sigma = Sigma_true)
##----------------------------------------------------------------##

##----------------------------------------------------------------##
## Create TMB objects
TMB_data = list(Y = Y)
TMB_params = list(logs_sd_RE=runif(n = d, min = 0, max = 2),
                   cov_RE= runif(n = ((d)*(d)-(d))/2, min = 0.1, max = 0.2))
##----------------------------------------------------------------##

##----------------------------------------------------------------##
## Run TMB
obj <- MakeADFun(data = TMB_data, parameters = TMB_params, DLL="tmb_MVN")
obj$hessian <- TRUE
opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max, trace=T)
# opt <- do.call("optim", obj)
# opt
# opt$hessian ## <-- FD hessian from optim
report = sdreport(obj)
##----------------------------------------------------------------##

##----------------------------------------------------------------##
## Get the estimated covariance matrices and check that they are PSD
## using the covariances as they are, and exponentiating logs_sd_RE for the diagonal
#

##' The way I populate the covariance matrix is as follows:
##' 1. I populate the off-diagonal entries as specified in 
##'    https://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html
##' 2. for the diagonal entries, first I exponentiate it (as I am estimating the
##'    log of the standard deviations to guarantee values > 0). I also take the
##'    square, as the inferred values when you use VECSCALE_t to scale the matrices
##'    are, apparently, the standard deviations (see TMB documentation above,
##'    although it seems contradictory: it says "Sigma has 1's on its diagonal
##'    To scale the variances we can use VECSCALE_t", but also " Set all standard
##'    deviations to 2.0")
cov_mat_est = fill_covariance_matrix(
       arg_d = d,
       arg_entries_var = exp(python_like_select_name(report$par.fixed, "logs_sd_RE"))**2,
       arg_entries_cov = python_like_select_name(report$par.fixed, "cov_RE"))

##' This covariance is not always PSD, however:

## this throws an error if the estimated sigma is not PSD
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est)

## another package to determine if a matrix is PSD (even if numerically isn't , exactly, I think!)
matrixcalc::is.positive.definite(cov_mat_est)

##' Here I started wondering if I was misinterpreting the way UNSTRUCTURED_CORR_t works,
##' and if I was populating it incorrectly

##' Alternative 1: without taking the square of the standard deviations for the diagonal
##' i.e. assuming that the parameter "logs_sd_RE" is actually the variances, and not the
##' standard deviations (which could well be, even though in the example on the TMB webpage
##' they call them "sds":
##' 
##' UNSTRUCTURED_CORR_t<Type> nll(theta); ## from TMB website, where theta is a vector of covariances
##' res = VECSCALE_t(nll,sds)(x);
##' 
##' in my case I do, equivalently,
##' 
##' UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_RE);
##' nll += VECSCALE_t(nll_mvn, exp(logs_sd_RE))(Y.row(i));
##' 
##' It would make sense that sds are actually variances, and not standard deviations,
##' as the example of TMB doesn't use a log-transformation

## The results of this new matrix are also not PSD
cov_mat_est_v2 <- NA
cov_mat_est_v2 = fill_covariance_matrix(arg_d = d,
  arg_entries_var = exp(python_like_select_name(report$par.fixed, "logs_sd_RE")),
  arg_entries_cov = python_like_select_name(report$par.fixed, "cov_RE"))
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_v2)

##' Alternative 2: throughout TMB they talk about "correlation matrix", not "covariance matrix".
##' First of all, I don't get how this could be a correlation matrix as the values are not in (-1,1)
##' 
##' In any case, in this second alternative I try to rescale the values cov_RE,
##' which I now take to be correlations, to get covariances
cov_mat_est_v3 = give_UNSTRUCTURED_CORR_t_matrix(python_like_select_name(report$par.fixed, "cov_RE"), dim_mat = d)
## for each element, get the covariance from the correlation, by multiplying by the standard deviations
.stdevs <- python_like_select_name(report$par.fixed, "logs_sd_RE")
for(i in 1:d){
  for(j in i:d){
    cov_mat_est_v3[i,j] = cov_mat_est_v3[j,i] =cov_mat_est_v3[i,j]*.stdevs[i]*.stdevs[j]
  }
}
## we don't need to add the variances in the diagonal because we already have them:
all(diag(cov_mat_est_v3) == .stdevs**2)

mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_v3)  ## not PSD

##' it doesn't look like the third alternative is the correct one.
##' my first option is what gives the estimated matrix most similar to the true matrix
par(mfrow=c(1,4))
image(Sigma_true, main='True')
image(cov_mat_est, main='logs_sd_RE are sds')
image(cov_mat_est_v2, main='logs_sd_RE are variances')
image(cov_mat_est_v3, main='cov_RE are corrs')

##' This pairs plot shows the same
my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}
pairs(cbind.data.frame(true=as.vector(Sigma_true), as.vector(cov_mat_est),
                        as.vector(cov_mat_est_v2), as.vector(cov_mat_est_v3)),
      lower.panel = my_line, upper.panel = my_line)

##' My question is: it looks like I populated the matrix correctly in cov_mat_est,
##' but the matrix is not PSD. What is going on? In quite some cases I do get estimates =
##' that are PSD

cov_mat_est_experimental <- cov_mat_est
diag(cov_mat_est_experimental) <- diag(cov_mat_est_experimental)*10

diag(Sigma_true)
Sigma_true
cov_mat_est

cov_mat_est_experimental <- cov_mat_est
diag(cov_mat_est_experimental) <- exp(python_like_select_name(report$par.fixed, "logs_sd_RE"))**4
cov_mat_est_experimental**2

mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_experimental)

cov_mat_est_x = fill_covariance_matrix(
  arg_d = d,
  arg_entries_var = (python_like_select_name(report$par.fixed, "logs_sd_RE"))**2,
  arg_entries_cov = python_like_select_name(report$par.fixed, "cov_RE"))
mvtnorm::rmvnorm(1, mean = rep(0,d), sigma = cov_mat_est_x)
