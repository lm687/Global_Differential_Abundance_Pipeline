## Assess how good the model is for determining differential abundance

rm(list = ls())

library(uuid)
source("2_inference/helper/helper_DA_stan.R")
source("1_create_ROO/roo_functions.R")
source("1_create_ROO/helper_1_create_ROO.R")


## (A) simulate from the model with different parameters
d = 5 # opt$Nk ## number of signatures
n = 100 #opt$Ns ## number of samples
beta_gamma_shape = 2.5  ##opt$hyperparam_shape ## shape parameter for the beta
sd_RE = 2 ## standard deviation for random effects
lambda = c(10, 10) ## overdispersion scalars. Lower value -> higher overdispersion
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
alphabar = softmax_mat( cbind(t(X_sim)%*%beta + t(Z_sim)%*%replicate(d-1, u, simplify = TRUE), 0) )
alpha = alphabar * lambdas

# image(t(alpha), main='Alphas')

## Create the counts
Nm = rpois(n*2, lambda = Nm_lambda) # number of mutations per sample \in N^{2*n}
W = matrix(NA, nrow = 2*n, ncol = d)
for(l in 1:(2*n)){
  W[l,] = HMP::Dirichlet.multinomial(Nrs = Nm[l], shape = alpha[l,])
}

image(t(W))

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

uuid = uuid::UUIDgenerate()
write.table("3_analysis/helper/table_simulation_params.txt", append = FALSE, x = cbind('d', 'n', 'beta_gamma_shape', 'sd_RE', 'Nm_lambda'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table("3_analysis/helper/table_simulation_params.txt", append = TRUE, x = cbind(d, n, beta_gamma_shape, sd_RE, Nm_lambda), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

saveRDS(list(objects_counts=objects_counts, d=d, n= n, beta_gamma_shape=beta_gamma_shape, sd_RE=sd_RE, lambda=lambda, Nm_lambda=Nm_lambda,
             X_sim = X_sim, beta = beta, Z_sim = Z_sim, u = u, lambdas = lambdas, alphabar = alphabar, alpha = alpha, Nm = Nm, W = W),
        file = paste0("../data/assessing_models_simulation/datasets/ ", uuid, ".RDS"))

## (B) We simulate two populations, and mix them with some proportions.


