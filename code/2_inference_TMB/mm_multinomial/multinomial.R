rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
compile("multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("multinomial"))
# set.seed(123)
n = 100
d = 4
true_theta = MCMCpack::rdirichlet(1, rep(1, d))
data <- list(Y = t(rmultinom(n, 300, true_theta)), n=n, d=d)
parameters <- list(theta_prime=rep(1, d-1))
obj <- MakeADFun(data, parameters, DLL="multinomial")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
# obj$he()    ## <-- Analytical hessian
sdreport(obj)

closure = function(v) v/sum(v)

theta_est = closure(exp(c( opt$par, 0)))

plot(true_theta, theta_est)

