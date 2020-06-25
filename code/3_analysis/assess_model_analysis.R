
# scp -r morril01@clust1-headnode.cri.camres.org:/Users/morril01/git_phd/out/DA_stan/DM_B_* /Users/morril01/Documents/PhD/CDA_in_Cancer/out/DA_stan/

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(rstan)
library(optparse)
library(ggplot2)
library(cowplot)
source("helper/helper_analyse_posteriors.R")

flder_inference = "../../data/assessing_models_simulation/inference_results/"
flder_datasets = "../../data/assessing_models_simulation/datasets/"
fles = list.files(flder_inference, full.names = TRUE)
fles = fles[grep("_DM.RDS", fles)]
names_uuid = strsplit(basename(fles), '_')[[1]][1]
names(names_uuid) = fles

posteriors = list()
for(f in fles){
  load(f)
  # cat(f, '\n')
  posteriors[[f]] = tryCatch(extract(fit))
}

## checking if there are any files in which the sampling didn't work
check_if_no_samples(posteriors, fles, flder_inference)

#pdf("~/Desktop/DM_all_by_it.pdf")
par(mfrow=c(1,length(alpha)))
posteriors = true_beta = in_cred_int_alpha = opt_params = list()
for(f in fles){
  load(f)
  dataset = readRDS(paste0(flder_datasets, names_uuid[[f]], '.RDS'))
  if(!exists("opt")){
    ## those runs were prior to having arguments for Nk, Ns, Nm_lambda
    next
  }else{
    posteriors[[f]] = tryCatch(rstan::extract(fit_stan))
    
    opt_params[[f]] = c(dataset$d, dataset$n, dataset$Nm_lambda)
    true_beta[[f]] = dataset$beta
    
    ## we select [,2,] to get only the slope, not the intercept (which by default is zero)
    in_cred_int_alpha[[f]] = any(sapply(1:(d-1), function(nk_it)
      get_bool_in_credible_interval(posterior = posteriors[[f]]$beta[,2,][,nk_it],
                                    true = 0))) ## is it differentially abundant?
    percentage_success_cred_int = 1 - (in_cred_int_alpha[[f]] - any(dataset$beta[2,] != 0)) ## if 1, success. if 0, failure
  }
}
dev.off()#

percentage_success_cred_int = sapply(in_cred_int_alpha, function(i) sum(i)/length(i))

par(mfrow=c(1,3))
# posteriors = true_params =  =  = list()
# plot(sapply(opt_params, function(i) i[1]),
#      percentage_success_cred_int)
# plot(sapply(opt_params, function(i) i[2]),
#      percentage_success_cred_int)
# plot(sapply(opt_params, function(i) i[3]),
#      percentage_success_cred_int)

df_multinomialrecovery = do.call('rbind', lapply(opt_params, rbind))
df_multinomialrecovery = cbind(df_multinomialrecovery, percentage_success_cred_int)
colnames(df_multinomialrecovery) = c('Nk', 'Ns', 'Nm_lambda', 'percentage_success_cred_int')
rownames(df_multinomialrecovery) = gsub(paste0(flder_inference, '/', collapse=''), '', rownames(df_multinomialrecovery))
rownames(df_multinomialrecovery) = gsub('.Rdata', '', rownames(df_multinomialrecovery))
dim(df_multinomialrecovery)

df_multinomialrecovery = as.data.frame(df_multinomialrecovery)
titl = ggdraw() +
  draw_label(
    "Accuracy in parameter recovery (percentage within credible interval)",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(titl,
          plot_grid(ggplot(as.data.frame(df_multinomialrecovery), aes(x=Nk, group=Nk, y=percentage_success_cred_int))+geom_violin()+geom_jitter()+ggtitle('Nk')+labs(y='% sigs correct alpha'),
                    #ggplot(as.data.frame(df_multinomialrecovery), aes(x=Nk, group=Nk, y=percentage_success_cred_int))+geom_boxplot()+ggtitle('Nk')+labs(y='% sigs correct alpha'),
                    ggplot(as.data.frame(df_multinomialrecovery), aes(x=Ns, group=Ns, y=percentage_success_cred_int))+geom_violin()+geom_jitter()+ggtitle('Ns')+labs(y='% sigs correct alpha'),
                    #ggplot(as.data.frame(df_multinomialrecovery), aes(x=Ns, group=Ns, y=log(percentage_success_cred_int)))+geom_boxplot()+ggtitle('Ns')+labs(y='% sigs correct alpha'),
                    ggplot(as.data.frame(df_multinomialrecovery), aes(x=Nm_lambda, group=Nm_lambda, y=log(percentage_success_cred_int)))+geom_violin()+geom_jitter()+ggtitle('Nm lambda')+labs(y='% sigs correct alpha'),
                    #ggplot(as.data.frame(df_multinomialrecovery), aes(x=Nm_lambda, group=Nm_lambda, y=percentage_success_cred_int))+geom_boxplot()+ggtitle('Nm lambda')+labs(y='% sigs correct alpha'),
                    #labels = c('A', 'B'),
                    label_size = 12, ncol = 3), ncol=1, rel_heights = c(0.1, 1))
ggsave("../../../../results/ProjectDA/stan_DA/recovery_params_DM2.pdf", width = 10, height = 4)

## Example of one run (the last one)
par(mfrow=c(1,3))
plot_posterior_and_true(posterior = posteriors[[f]]$alpha, true = alpha,
                        main=gsub('.Rdata', '', gsub(paste0(flder_inference, '/'), '', f)),
                        default_par=FALSE)

colMeans(posteriors[[f]]$alpha)

fit

## recovery of beta
posteriors[[f]]

