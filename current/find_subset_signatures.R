rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/")
library(TMB)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/fullRE_ME_dirichletmultinomial"))
TMB::compile("mm_multinomial/fullRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_dirichletmultinomial_single_lambda"))
TMB::compile("mm_multinomial/fullRE_dirichletmultinomial_single_lambda2.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_dirichletmultinomial_single_lambda2"))
TMB::compile("mm_multinomial/diagRE_dirichletmultinomial_single_lambda.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_dirichletmultinomial_single_lambda"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/FE_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/FE_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/diagRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/diagRE_ME_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda2.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomialsinglelambda2"))
TMB::compile("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov.cpp", "-std=gnu++17")
dyn.load(dynlib("../../current/Dirichlet_Multinomial_Dom/code/TMB_models/sparseRE_ME_dirichletmultinomial_singlecov"))


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

res_FE_DM <- wrapper_run_TMB(object = Kidney_RCC_clearcell, model = "FE_DM")
res_FE_DM ## good convergence

res_diagME_DM_nlminb <- wrapper_run_TMB_debug(Kidney_RCC_clearcell, iter.max = 250, model = "diagRE_DM", return_report=T)
res_diagME_DM_nlminb ## relative convergence
res_diagME_DM <- wrapper_run_TMB(object = Kidney_RCC_clearcell, model = "diagRE_DM") ## seems to be much slower
wald_TMB_wrapper(res_diagME_DM)
wald_TMB_wrapper(res_diagME_DM_nlminb)

plot(res_diagME_DM$par.fixed, res_diagME_DM_nlminb$par, col=factor(is_slope(names(res_diagME_DM_nlminb$par))))
abline(coef = c(0,1), lty='dashed')

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

singlelambda_res <- wrapper_run_TMB_debug(object = Kidney_RCC_clearcell, iter.max = 500, init_log_lambda = NA, return_report = T,
                                          model = "fullREDMsinglelambda")
doublelambda_res <- wrapper_run_TMB_debug(object = Kidney_RCC_clearcell, iter.max = 500, init_log_lambda = NA, return_report = T,
                                          model = "fullRE_DM")

singlelambda_res <- wrapper_run_TMB_debug(object = give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS12', 'SBS41', 'SBS2')),
                                          iter.max = 500, init_log_lambda = NA, return_report = T,
                                          model = "fullREDMsinglelambda")
singlelambda_res

##' give the indices of the covariance matrix (non-redudant lower diagonal entries) which we are interested in inferring
##' (i.e. those for non-zero exposures)
##' we are always going to sort the matrix by increasing order of total exposure, so that the baseline signature is one which
##' is clearly nonzero



#----------------------------------------------------------------------#

enough_samples = readLines("~/Desktop/CT_sufficient_samples.txt")
df_all_samples = data.frame(do.call('rbind', lapply(enough_samples, function(i) rbind(c(i, 'signatures'), c(i, "nucleotidesubstitution1")))))
nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t", comment.char = "#", fill = F)
library(parallel)

mclapply(1:nrow(df_all_samples),
         function(idx){
           i = df_all_samples[idx,]
           x = wrapper_run_TMB_debug(object = load_PCAWG(ct = i[1,1], typedata = i[1,2]),
                               model = "diagREDMsinglelambda", return_report=T)
           saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "diagRE_DMSL_",
                                           paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
         })
ct = "Bone-Osteosarc"
.x <- wrapper_run_TMB(object = load_PCAWG(ct = ct, typedata = 'signatures'), model = "diagREDMsinglelambda")
saveRDS(object = .x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "diagRE_DMSL_",
                                paste0(ct, "_signatures", collapse = "_"), ".RDS"))

#----------------------------------------------------------------------#

subset_sigs <- read.table("../../current/subset_sigs.txt", stringsAsFactors = F, fill = T)
.pcawg0 <- load_PCAWG(ct = ct, typedata = "signatures")
# .sigs_rm <- setdiff(colnames(.pcawg$Y), strsplit(subset_sigs[subset_sigs$V2 == ct,'V3'], ',')[[1]])
.sigs_rm <- 'SBS33' #colnames(.pcawg0$Y)[-c(1:4)]
.pcawg <- give_subset_sigs_TMBobj(.pcawg0, sigs_to_remove = .sigs_rm)
res_subset <- wrapper_run_TMB(object = sort_columns_TMB(.pcawg),
                                     model = "fullREDMsinglelambda2", use_nlminb = T)
res_subset
# wald_TMB_wrapper(res_sparse)
saveRDS(res_subset, paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/subsetsorted_RE_DMSL_", ct, "_signatures.RDS"))

# mclapply(which(df_all_samples$X2 == "signatures"),
#          function(idx){
#            i = df_all_samples[idx,]
#            typedata = i[1,2]
#            obj_subset <- give_subset_sigs_TMBobj(load_PCAWG(ct = i[1,1], typedata = i[1,2]),
#                                                  sigs_to_remove = unique(nonexogenous$V1))
#            if(dim(obj_subset$Y)[2] <  dim(load_PCAWG(ct = i[1,1], typedata = i[1,2])$Y)[2]){
#              x = wrapper_run_TMB_debug(object = obj_subset, 
#                                        model = "fullRE_DM", return_report=T)
#              x
#              saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fulLRE_nonexo_DM_",
#                                              paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
#            }
#          })

#----------------------------------------------------------------------#

# give_subset_sigs_roo = function(roo, subset_sigs, slot_name){
#   #subset_sigs=c('SBS1', 'SBS2', 'SBS5', 'SBS13', 'SBS18', 'SBS40')
#   # roo=readRDS(paste0("../../data/", "roo/", ct, '_', "signatures", "_ROO.RDS" ))
#   roo
# }

enough_samples = readLines("~/Desktop/CT_sufficient_samples.txt")
subset_sigs_sparse_cov_idx <- read.table("../../current/subset_sigs_sparse_cov_idx.txt", stringsAsFactors = F, fill = T)
# ct <- enough_samples[29]
# obj_ct = load_PCAWG(ct, typedata="signatures", path_to_data = "../../data/")
# idx_cov_to_fill_read <- as.integer(strsplit(subset_sigs_sparse_cov_idx[subset_sigs_sparse_cov_idx$V1 == ct,2], ",")[[1]])

# idx_cov_to_fill_read = idx_cov_to_fill_read[-1]
# idx_cov_to_fill_read = idx_cov_to_fill_read[1:2]
# res_sparse <- wrapper_run_TMB_debug(object = obj_ct, iter.max = 500, init_log_lambda = 3, return_report = T,
                                    # model = "fullREDMsinglelambda") ## idx_cov_to_fill must be zero-indexed
# res_sparse
# xxxx <- wrapper_run_TMB(object = sort_columns_TMB(obj_ct), model = "fullREDMsinglelambda2", use_nlminb = T, smart_init_vals = T)
res_sparse2 <- wrapper_run_TMB_debug(object = sort_columns_TMB(obj_ct),
                                     iter.max = 500, init_log_lambda = 3, return_report = T,
                                    model = "sparseRE_DMSL2",
                                    idx_cov_to_fill=idx_cov_to_fill_read-1) ## idx_cov_to_fill must be zero-indexed
res_sparse2
# wald_TMB_wrapper(res_sparse)
# saveRDS(res_sparse, paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL_", ct, "_signatures.RDS"))
saveRDS(res_sparse2, paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/sparseRE_DMSL2_", ct, "_signatures.RDS"))
#----------------------------------------------------------------------#

res_DL2 <- wrapper_run_TMB(object = obj_ct,
                                    model = "fullREDMsinglelambda2")
res_DL2

mclapply(which(df_all_samples$X2 == "signatures"),
         function(idx){
           i = df_all_samples[idx,]
           typedata = i[1,2]
           obj_subset <- give_subset_sigs_TMBobj(load_PCAWG(ct = i[1,1], typedata = i[1,2]),
                            sigs_to_remove = unique(nonexogenous$V1))
           if(dim(obj_subset$Y)[2] <  dim(load_PCAWG(ct = i[1,1], typedata = i[1,2])$Y)[2]){
             x = wrapper_run_TMB_debug(object = obj_subset, 
                                       model = "fullRE_DM", return_report=T)
             x
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fulLRE_nonexo_DM_",
                                             paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
             # idx_cov_to_fill_read <- as.integer(strsplit(subset_sigs_sparse_cov_idx[subset_sigs_sparse_cov_idx$V1 == ct,2], ",")[[1]])
             # x = wrapper_run_TMB_debug(object = obj_subset, iter.max = 500, init_log_lambda = 3, return_report = T,
             #                       model = "sparseRE_DM",
             #                       idx_cov_to_fill=idx_cov_to_fill_read-1)
             # saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "sparseRE_nonexo_DM_",
             #                                 paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
           }
         })

mclapply(which(df_all_samples$X2 == "signatures"),
         function(idx){
           i = df_all_samples[idx,]
           typedata = i[1,2]
           obj_subset <- give_subset_sigs_TMBobj(load_PCAWG(ct = i[1,1], typedata = i[1,2]),
                                                 sigs_to_remove = unique(nonexogenous$V1))
           if(dim(obj_subset$Y)[2] <  dim(load_PCAWG(ct = i[1,1], typedata = i[1,2])$Y)[2]){
             x = wrapper_run_TMB_debug(object = obj_subset, 
                                       model = "fullREDMsinglelambda", return_report=T)
             x
             saveRDS(object = x, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fullRE_nonexo_DMSL_",
                                             paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
           }
         })

subset_sigs_sparse_cov_idx_nonexo <- read.table("../../current/subset_sigs_sparse_cov_idx_nonexo.txt", stringsAsFactors = F, fill = T)

mclapply(which(df_all_samples$X2 == "signatures"),
         function(idx){
           i = df_all_samples[idx,]
           typedata = i[1,2]
           ct = i[1,1]
           obj_subset <- give_subset_sigs_TMBobj(load_PCAWG(ct = i[1,1], typedata = i[1,2]),
                                                 sigs_to_remove = unique(nonexogenous$V1))
           idx_cov_to_fill_read <- as.integer(strsplit(subset_sigs_sparse_cov_idx_nonexo[subset_sigs_sparse_cov_idx_nonexo$V1 == ct,2], ",")[[1]])
           
           idx_cov_to_fill_read
           if(dim(obj_subset$Y)[2] <  dim(load_PCAWG(ct = i[1,1], typedata = i[1,2])$Y)[2]){
             .res_sparse <- wrapper_run_TMB_debug(object = sort_columns_TMB(obj_subset),
                                                  iter.max = 500, init_log_lambda = 3, return_report = T,
                                                  model = "sparseRE_DMSL2",
                                                  idx_cov_to_fill=idx_cov_to_fill_read-1) ## idx_cov_to_fill must be zero-indexed
             .res_sparse
             saveRDS(object = .res_sparse, file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "sparseRE_nonexo_DMSL_",
                                             paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
           }
         })
obj_subset2 <- give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1))
sort(colSums(obj_subset2$Y))
createBarplot(normalise_rw(makenames_row(obj_subset2$Y)))


aa <- readRDS(file=paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fulLRE_nonexo_DM_",
                                paste0(df_all_samples[idx,], collapse = "_"), ".RDS"))
wald_TMB_wrapper(aa)

fles_nonexo <- apply(df_all_samples[df_all_samples$X2 == "signatures",], 1, function(i){
  paste0("../../data/pcawg_robjects_cache/tmb_results/nlminb/", "fulLRE_nonexo_DM_",
         paste0(i, collapse = "_"), ".RDS")
})

t(sapply(fles_nonexo, function(j){
  aa <- try(readRDS(file=j))
  if(typeof(aa) == "character"){
    c(basename(j), NA)
  }else{
    c(j, wald_TMB_wrapper(aa))
  }
}))

makenames_row = function(a){
  rownames(a) = make.names(rownames(a), unique = T)
  a
}


createBarplot(normalise_rw(makenames_row(load_PCAWG(ct = ct, typedata = 'signatures')$Y)))
createBarplot(normalise_rw(makenames_row(obj_subset$Y)))
paste0(colnames(obj_subset$Y), collapse = ',')

ct="CNS-PiloAstro"
ct="CNS-Medullo"
ct <- "Prost-AdenoCA"
ct <- 'Lymph-BNHL'
ct <- 'LUNG-SCC'
ct <- 'Kidney-ChRCC'
obj_subset2 <- give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1))
sort(colSums(obj_subset2$Y))
createBarplot(normalise_rw(makenames_row(obj_subset2$Y)))

pdf("../../../CDA_in_Cancer/text/presentations/20210421_WLS_Morrill/figures/prostateadenoca_nonexo.pdf", width=9, height = 4)
obj_subset2 <- give_subset_sigs_TMBobj(load_PCAWG(ct = "Prost-AdenoCA", typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1))
do.call('grid.arrange', list(grobs=lapply(split.data.frame(normalise_rw(makenames_row(obj_subset2$Y)), f = rep(c(1,2), each=(nrow(obj_subset2$Y)/2))),
       function(i) (createBarplot(i, remove_labels = T))), ncol=2))
dev.off()

pdf("../../../CDA_in_Cancer/text/presentations/20210421_WLS_Morrill/figures/prostateadenoca_nonexo.pdf", width=9, height = 4)
obj_subset2 <- give_subset_sigs_TMBobj(load_PCAWG(ct = "CNS-PiloAstro", typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1))
do.call('grid.arrange', list(grobs=lapply(split.data.frame(normalise_rw(makenames_row(obj_subset2$Y)), f = rep(c(1,2), each=(nrow(obj_subset2$Y)/2))),
                                          function(i) (createBarplot(i, remove_labels = T))), ncol=2))
dev.off()


num_cats_nonexo <- sapply(unique(df_all_samples$X1), function(ct){
  .xx <- try(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                       sigs_to_remove = unique(nonexogenous$V1)))
  try(ncol(.xx$Y))
})
num_cats_nonexo
