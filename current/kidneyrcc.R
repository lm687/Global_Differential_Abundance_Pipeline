##' case of kidney cancer rcc
##' object=give_subset_sigs_TMBobj(Kidney_RCC_clearcell, sigs_to_remove = c('SBS2', 'SBS12', 'SBS41')) #"relative convergence (4)"

## sort columns
sorted_cols = names(sort(colSums(object$Y), decreasing = F))
# res <- wrapper_run_TMB_debug(object, Kidney_RCC_clearcell)

covmat = res$cov.fixed
rownames(covmat)[grep('beta', rownames(covmat))] = paste0('beta_', rep(sorted_cols[-length(sorted_cols)], each=2), '_',
                                                          rep(c('intercept', 'slope'), length(sorted_cols)-1))
colnames(covmat)[grep('beta', colnames(covmat))] = paste0('beta_', rep(sorted_cols[-length(sorted_cols)], each=2), '_',
                                                          rep(c('intercept', 'slope'), length(sorted_cols)-1))

cov_names <- paste0('cov_', outer(sorted_cols[-length(sorted_cols)], sorted_cols[-length(sorted_cols)], FUN = Vectorize(function(i, j){paste0(i, '-', j)}))[lower.tri(outer(sorted_cols[-length(sorted_cols)], sorted_cols[-length(sorted_cols)], FUN = Vectorize(function(i, j){paste0(i, '-', j)})))])
rownames(covmat)[grep('cov_par_RE', rownames(covmat))] = cov_names
colnames(covmat)[grep('cov_par_RE', colnames(covmat))] = cov_names

rownames(covmat)[grep('logs_sd_RE ', rownames(covmat))] = paste0('logs_', rep(sorted_cols[-length(sorted_cols)]))
colnames(covmat)[grep('logs_sd_RE ', colnames(covmat))] = paste0('logs_', rep(sorted_cols[-length(sorted_cols)]))

# pdf("../../current/kidney/covmat.pdf", height = 8, width = 8)                    
pheatmap::pheatmap(covmat)
# dev.off()

problematic_cov = python_like_select_rownames(summary(res), "cov_par_RE")
problematic_cov= cbind.data.frame(problematic_cov,
                                  name=paste0('cov_', outer(sorted_cols[-length(sorted_cols)], sorted_cols[-length(sorted_cols)], FUN = Vectorize(function(i, j){paste0(i, '-', j)}))[lower.tri(outer(sorted_cols[-length(sorted_cols)], sorted_cols[-length(sorted_cols)], FUN = Vectorize(function(i, j){paste0(i, '-', j)})))]))

# cov_par_RE.12  45.19582846        NaN  cov_SBS29-SBS1
# cov_par_RE.13  13.32040653        NaN   cov_SBS5-SBS1
# cov_par_RE.14  27.70145184        NaN  cov_SBS22-SBS1
# cov_par_RE.15   7.35778203        NaN  cov_SBS29-SBS6
# cov_par_RE.16  73.25063914        NaN   cov_SBS5-SBS6
# cov_par_RE.17  -8.28140191        NaN  cov_SBS22-SBS6
# cov_par_RE.18 -41.54147050        NaN  cov_SBS5-SBS29
# cov_par_RE.19  68.04243028        NaN cov_SBS22-SBS29

problematic_cov[is.na(problematic_cov$`Std. Error`),]
# Estimate Std. Error            name
# cov_par_RE.6 -2.799380        NaN cov_SBS22-SBS13
# cov_par_RE.8 -3.426543        NaN   cov_SBS1-SBS2
# cov_par_RE.9 -5.045325        NaN   cov_SBS6-SBS2

# Estimate Std. Error           name
# cov_par_RE.14 0.06958914        NaN cov_SBS6-SBS41
# cov_par_RE.25 0.20463902        NaN cov_SBS5-SBS29

# cov_par_RE.6 -2.227599        NaN cov_SBS22-SBS13
# cov_par_RE.7  2.195222        NaN  cov_SBS41-SBS2
# cov_par_RE.8 -2.594447        NaN   cov_SBS1-SBS2
# cov_par_RE.9 -3.724536        NaN   cov_SBS6-SBS2