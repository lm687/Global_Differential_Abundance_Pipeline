#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
# set.seed(1245)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
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

rsq = function (x, y) cor(x, y) ^ 2
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
compile("ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("ME_multinomial"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
d = 4 # opt$Nk ## number of signatures
n = 10 #opt$Ns ## number of samples
beta_gamma_shape = 2  ##opt$hyperparam_shape ## shape parameter for the beta
sd_RE = 0.3 ## standard deviation for random effects
lambda = c(200, 200) ## overdispersion scalars. Lower value -> higher overdispersion
Nm_lambda = 300 ## opt$Nm_lambda ## lambda parameter for number of mutations per sample (i.e. a sample in a group)

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
theta = softmax( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )

W = t(apply(theta, 1, rmultinom, n = 1, size = Nm_lambda))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
data <- list(Y = W, n=n*2, d=d, x=t(X_sim), z = t(Z_sim), num_individuals=n)
# dev.off(); image(t(data$Y))

parameters <- list(
  beta = array(c(rep(runif(1, min = -4, max = 4), (d-1)),
                 rep(runif(1, min = -4, max = 4), (d-1))),
               dim = c(2,d-1)),
  u = matrix(rep(0, n)),
  logSigma_RE=0
)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
obj <- MakeADFun(data, parameters, DLL="ME_multinomial")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
## Comparing true estimate values
par(mfrow=c(1,2))
plot(beta, opt$par[grep(names(opt$par), pattern = 'beta')])
abline(a = c(0,1))
plot(u, opt$par[grep(names(opt$par), pattern = 'u')])
rsq(x = as.vector(u),  y = opt$par[grep('u', names(opt$par))])
abline(a = c(0,1))
rsq(x = as.vector(beta),  y = opt$par[grep('beta', names(opt$par))])
#-------------------------------------------------------------------------------------------------#

