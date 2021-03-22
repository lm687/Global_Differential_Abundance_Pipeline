rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/")
library(TMB)
library(ggplot2)
library(dplyr)
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/fullRE_ME_dirichletmultinomial"))

# TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
# dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial"))

unique(sapply(list.files("../../data/roo/"), function(i) strsplit(i, '_')[[1]][1]))
Kidney_RCC_clearcell = load_PCAWG("Kidney-RCC.clearcell", typedata="signatures", path_to_data = "../../data/")
CNS_Medullo = load_PCAWG("CNS-Medullo", typedata="signatures", path_to_data = "../../data/")


Kidney_RCC_clearcell$Y[,c('SBS12', 'SBS41')]
xxx = give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS12', 'SBS41'))
sort(colSums(Kidney_RCC_clearcell$Y))
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS12', 'SBS41')) #Error in if (m < 0) { : missing value where TRUE/FALSE needed
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS2', 'SBS41')) #"iteration limit reached without convergence (10)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS2', 'SBS12')) #"singular convergence (7)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS2', 'SBS12', 'SBS41')) #"relative convergence (4)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS2', 'SBS12', 'SBS41', 'SBS1')) #"relative convergence (4)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS13', 'SBS2', 'SBS12', 'SBS41', 'SBS1')) #"relative convergence (4)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS13', 'SBS2', 'SBS12', 'SBS41', 'SBS1', 'SBS6')) #  "relative convergence (4)"
# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('dummy'))
# 
# res <- wrapper_run_TMB_debug(object, Kidney_RCC_clearcell, iter.max = 250)
# 
# 
# cat(paste0(colnames(Kidney_RCC_clearcell$Y)[!(colnames(Kidney_RCC_clearcell$Y) %in% c('dummy'))], collapse=', '))
# 
# all_combinations = do.call('c', lapply(1:(ncol(Kidney_RCC_clearcell$Y)-2), function(k){
#   .x = combn(colnames(Kidney_RCC_clearcell$Y), k)
#   lapply(seq_len(ncol(.x)), function(i) .x[,i])
#   }))
# all_combinations
# 
# idx_combs = 1:80
# res_TMB = lapply(all_combinations[idx_combs], function(comb_exclude) wrapper_run_TMB_debug(give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = comb_exclude), Kidney_RCC_clearcell) )
# cat((apply(cbind(sapply(all_combinations[idx_combs], function(j) paste0(colnames(Kidney_RCC_clearcell$Y)[!(colnames(Kidney_RCC_clearcell$Y) %in% j)], collapse=', ')),
#                     sapply(res_TMB, `[`, "message")), 1, paste0, collapse=': ')), sep = '\n')

# opt$hessian ## <-- FD hessian from optim
# aa = sdreport(obj)
# aa

# do.call('grid.arrange', lapply(split_Y(give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS23', 'SBS39', 'SBS40'))$Y),
#                                function(i) print(createBarplot(normalise_rw(i)))))
# do.call('grid.arrange', lapply(split_Y(give_subset_sigs_TMBobj(Kidney_RCC_clearcell,
#                                                                sigs_to_remove = colnames(Kidney_RCC_clearcell$Y)[!(colnames(Kidney_RCC_clearcell$Y) %in% c('SBS23', 'SBS39', 'SBS40'))])$Y),
#                                function(i) print(createBarplot(normalise_rw(i)))))

TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/FE_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/FE_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/diagRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/diagRE_ME_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov"))

res_FE_DM <- wrapper_run_TMB(object = Kidney_RCC_clearcell, model = "FE_DM")
res_FE_DM ## good convergence

res_diagME_DM <- wrapper_run_TMB_debug(Kidney_RCC_clearcell, iter.max = 250, model = "diagRE_DM")
res_diagME_DM ## relative conervergence

# object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS12'))
res <- wrapper_run_TMB_debug(Kidney_RCC_clearcell, iter.max = 500, init_log_lambda = 3, return_report = T)
res

plot(res$gradient.fixed, type='h', col=factor(is.na(python_like_select_rownames(summary(res), c('beta|cov_par_RE|log_lambda'))[,2])))

res_sparse <- wrapper_run_TMB_debug(object = Kidney_RCC_clearcell, iter.max = 500, init_log_lambda = 3, return_report = T, model = "sparseRE_DM",
                                    idx_cov_to_fill=c(2L, 4L, 3L, 8L, 7L, 15L)-1) ## idx_cov_to_fill must be zero-indexed
res_sparse
saveRDS(res_sparse, "../../data/pcawg_robjects_cache/tmb_results/sparseRE_DM_Kidney-RCC.clearcell_signatures.RDS")

res_sparse_CNSmed <- wrapper_run_TMB_debug(object = CNS_Medullo, iter.max = 500, init_log_lambda = 3, return_report = T, model = "sparseRE_DM",
                                    idx_cov_to_fill=c(23L, 22L, 21L, 24L, 12L, 11L, 14L, 6L, 9L, 4L)-1) ## idx_cov_to_fill must be zero-indexed
res_sparse_CNSmed

object = Kidney_RCC_clearcell; iter.max = 500; init_log_lambda = 3; return_report = T; model = "sparseRE_DM"; idx_cov_to_fill=c(2, 3, 6)

##' give the indices of the covariance matrix (non-redudant lower diagonal entries) which we are interested in inferring
##' (i.e. those for non-zero exposures)
##' we are always going to sort the matrix by increasing order of total exposure, so that the baseline signature is one which
##' is clearly nonzero
