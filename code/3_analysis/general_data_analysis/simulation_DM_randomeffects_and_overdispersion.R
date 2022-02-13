## simulated data with different implementations of the overdispersion parameter
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

set.seed(1234)

library(compositions)
library(HMP)
source("../../2_inference_TMB/helper_TMB.R")

## simulate some correlation structure. There are 3 signatures

a <- runif(100); b <- a+runif(100, min = -0.4, max=0.2); sigma <- cov(cbind(a,b))
sigma

beta <- rbind(c(0.5, 0.6), c(0.3, -1))

## simulate log-ratios
n <- 1000

## sample random intercepts
u <- mvtnorm::rmvnorm(n, mean=c(0,0), sigma = sigma)

x <- cbind(1, c(rep(0, n),rep(1, n)))
z <- give_z_matrix(n*2)
logR <- x %*% beta + z %*% u
probs <- softmax(cbind(logR, 0))

plot(compositions::acomp(probs), col=factor(x[,2]))

Nm <- rpois(n*2, lambda = 90)
Nm_b <- rpois(n*2, lambda = 25)
hist(Nm)

## add shared overdispersion
lambda <- c(200,200)
lambda <- c(20,20)
alpha <- rep(lambda, each=n)*probs

## add two overdispersion parameters
lambda_b <- c(9,6)
alpha_OVERDISPb <- rep(lambda_b, each=n)*probs

lambda_c <- c(200,6)
alpha_OVERDISPc <- rep(lambda_c, each=n)*probs


## get DM counts
give_DM <- function(alpha, Nm){
  t(sapply(1:nrow(alpha), function(i){
  HMP::Dirichlet.multinomial(Nrs = Nm[i], shape = alpha[i,] )
}))
}


## modify sigma
sigma_SIGMAb <- sigma *9
sigma_SIGMAc <- sigma *0.4
beta_b <- beta/9
u_SIGMAb <- mvtnorm::rmvnorm(n, mean=c(0,0), sigma = sigma_SIGMAb)
u_SIGMAc <- mvtnorm::rmvnorm(n, mean=c(0,0), sigma = sigma_SIGMAc)
logR_SIGMAb <- x %*% beta + z %*% u_SIGMAb
logR_SIGMAb_BETAb <- x %*% beta_b + z %*% u_SIGMAb
logR_SIGMAc_BETAb <- x %*% beta_b + z %*% u_SIGMAc
logR_BETAb <- x %*% beta_b + z %*% u
probs_SIGMAb <- softmax(cbind(logR_SIGMAb, 0))
probs_SIGMAb_BETAb <- softmax(cbind(logR_SIGMAb_BETAb, 0))
probs_SIGMAc_BETAb <- softmax(cbind(logR_SIGMAc_BETAb, 0))
probs_BETAb <- softmax(cbind(logR_BETAb, 0))
alpha_SIGMAb <- rep(lambda, each=n)*probs_SIGMAb
alpha_SIGMAb_BETAb <- rep(lambda, each=n)*probs_SIGMAb_BETAb
alpha_SIGMAc_BETAb <- rep(lambda, each=n)*probs_SIGMAc_BETAb
alpha_BETAb <- rep(lambda, each=n)*probs_BETAb

W <- give_DM(alpha, Nm)
W_NMb <- give_DM(alpha, Nm_b)
W_OVERDISPb <-  give_DM(alpha_OVERDISPb, Nm)
W_OVERDISPc <-  give_DM(alpha_OVERDISPc, Nm)
W_SIGMAb <-  give_DM(alpha_SIGMAb, Nm)
W_SIGMAb_BETAb <-  give_DM(alpha_SIGMAb_BETAb, Nm)
W_SIGMAc_BETAb <-  give_DM(alpha_SIGMAc_BETAb, Nm)
W_BETAb <-  give_DM(alpha_BETAb, Nm)


## re-normalise counts

# Wnorm <- normalise_rw(W)
# Wnorm_OVERDISPb <- normalise_rw(W_OVERDISPb)


## plot

par(mfrow=c(4,4), mar=c(0,0,0,0))
pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_OVERDISPb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_OVERDISPb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_OVERDISPc.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_OVERDISPc), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_SIGMAb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_SIGMAb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_NMb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_NMb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_SIGMAb_BETAb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_SIGMAb_BETAb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_BETAb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_BETAb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_SIGMAc_BETAb.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_SIGMAc_BETAb), plot_ternary, legend_on=F, plot_points=F)
dev.off()

##' things that can vary
##' - sigma for the random effects
##' - ratio of beta and u (which, supposing they are constant, should give the same log-ratios) (!!!!!)
##' - Nm: number of mutations, which also determine the variance in the normalised count exposures (but which is fixed in the model)
##' - overdispersion parameter, and whether there are two

## compensating lower overdispersion by having a higher variance for the random effects
## initial situation: the one we have above
## second situation, with lambda*, Sigma, where lambda* >> lambda (i.e. less overdispersion)
## third situation, with lambda*, Sigma*, in which the (low) overdispersion of lambda* is compensated by larger RE, i,e.
## higher variance in Sigma

lambda_d <- c(6000,6000)
alpha_OVERDISPd <- rep(lambda_d, each=n)*probs
W_OVERDISPd <-  give_DM(alpha_OVERDISPd, Nm)

pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_OVERDISPd.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_OVERDISPd), plot_ternary, legend_on=F, plot_points=F)
dev.off()

# sigma_SIGMAd <- sigma
# diag(sigma_SIGMAd) <- diag(sigma_SIGMAd)*5
sigma_SIGMAd <- sigma*7
u_SIGMAd <- mvtnorm::rmvnorm(n, mean=c(0,0), sigma = sigma_SIGMAd)
logR_SIGMAd <- x %*% beta + z %*% u_SIGMAd
probs_SIGMAd <- softmax(cbind(logR_SIGMAd, 0))
alpha_SIGMAd_OVERDISPd <- rep(lambda_d, each=n)*probs_SIGMAd
W_SIGMAd_OVERDISPd <-  give_DM(alpha_SIGMAd_OVERDISPd, Nm)


pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_W_SIGMAd_OVERDISPd.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(W_SIGMAd_OVERDISPd), plot_ternary, legend_on=F, plot_points=F)
dev.off()

<<<<<<< HEAD
=======

## why we need two betas
beta[2,] <- c(0,0)
logR <- x %*% beta + z %*% u
probs <- softmax(cbind(logR, 0))
alpha_OVERDISPc <- give_DM(rep(c(800,30), each=n)*probs, Nm = rep(90, n*2))
pdf(file = "../../../results/models_explanatory/DM_parameters_and_identifiability_why_differential_precision.pdf", heigh=3)
par(mfrow=c(1,2), mar=c(0,0,0,0))
sapply(split_matrix_in_half(alpha_OVERDISPc), plot_ternary, legend_on=F, plot_points=F)
dev.off()

>>>>>>> b7516544d6581da5bf0a960e309788c1fba6dff6
