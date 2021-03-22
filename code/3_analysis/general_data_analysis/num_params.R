rm(list = ls())
setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
require(reshape2)
require(ggplot2)
require(ggrepel)
folder_objects="../data/roo/"
folder_robjs = "../data/robjects_cache/tmb_results/"
#source("mm_multinomial/helper_functions.R")
source("2_inference_TMB/helper_TMB.R")

#----------------------------------------------------------------------------------------#
## Read all the ROO files which contain exposures in two groups
fles_in = list.files(folder_objects, full.names=TRUE)
roo_files = sapply(fles_in, readRDS)
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
results_TMB_M = lapply( python_like_select(list.files(folder_robjs), "^M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_M) = sapply(python_like_select(list.files(folder_robjs), "^M_"), clean_name)
results_TMB_DM = lapply( python_like_select(list.files(folder_robjs), "^DM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_DM) = sapply(python_like_select(list.files(folder_robjs), "^DM_"), clean_name)
results_TMB_DM_dep = lapply( python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"),
                             function(i) readRDS(paste0("../../data/robjects_cache/tmb_results_dep/", i)))
names(results_TMB_DM_dep) = sapply(python_like_select(list.files("../../data/robjects_cache/tmb_results_dep/"), "^DM_"), clean_name)
results_TMB_LNM = lapply( python_like_select(list.files(folder_robjs), "^LNM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_LNM) = sapply(python_like_select(list.files(folder_robjs), "^LNM_"), clean_name)
results_TMB_fullRE_M = lapply( python_like_select(list.files(folder_robjs), "^fullRE_M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_M) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_M_"), clean_name_fullRE)
results_TMB_fullRE_DM = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DM_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DM) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DM_"), clean_name_fullRE)
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
M = function(d){d-1}
M_FE = function(d){2*(d-1)}
M_ME = function(d){1/2*d**2 + 1/2*d-1}
DM = function(d){d}
DM_FE = function(d){2*d}
DM_ME = function(d){1/2*d**2 + 1/2*d+1}
LNM = function(d){1/2*(d**2) + 1/2*d - 1}
LNM_FE = function(d){1/2*d**2 + 5/2*d - 3}
LNM_ME = function(d){d**2 + 4*d - 2}
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
sq = 1:20
x = melt(cbind(M=sapply(sq, M),
      M_FE=sapply(sq, M_FE),
      M_ME=sapply(sq, M_ME),
      DM=sapply(sq, DM),
      DM_FE=sapply(sq, DM_FE),
      DM_ME=sapply(sq, DM_ME),
      LNM=sapply(sq, LNM),
      LNM_FE=sapply(sq, LNM_FE),
      LNM_ME=sapply(sq, LNM_ME)))
x[,'Model'] = sapply(as.character(x$Var2), function(i) strsplit(i, '_')[[1]][1])
x[,'Effects'] = sapply(as.character(x$Var2), function(i) strsplit(i, '_')[[1]][2])
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
get_sigs = function(i){
  if(is.na(i) | typeof(i) == "character"){
    NA
  }else{
    if(sum(sapply(i@count_matrices_active, length)) == 0){
      ncol(i@count_matrices_all[[1]])
    }else{
      ncol(i@count_matrices_active[[1]])
    }
  }
}
sapply(roo_files, get_sigs)
#----------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------#
names_sigs = basename(sort(names(roo_files[grepl('signatures', names(roo_files))])))
names_sigs_short = gsub("_", "", gsub("_ROO.RDS", "", names_sigs))
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# Single intercept
df_roo_sigs = data.frame(name=sapply(names_sigs,
                                     function(i) strsplit(i, '_')[[1]][1]),
                         num_samples=as.numeric(sapply(roo_files[sort(names(roo_files[grepl('signatures', names(roo_files))]))],
                                                       function(i) try(nrow(i@count_matrices_all[[1]])))),
                         num_sigs=sapply(roo_files[sort(names(roo_files[grepl('signatures', names(roo_files))]))], get_sigs),
                         convergence_M = sapply(results_TMB_M[python_like_select(names(results_TMB_M), 'signatures')[match(names_sigs_short,
                                                                                                                           names(python_like_select_name(results_TMB_M, 'signatures')))]], give_summary_per_sample),
                         convergence_DM = sapply(results_TMB_DM[python_like_select(names(results_TMB_DM), 'signatures')[match(names_sigs_short,
                                                                                                                              names(python_like_select_name(results_TMB_DM, 'signatures')))]], give_summary_per_sample),
                         convergence_LNM = sapply(results_TMB_LNM[python_like_select(names(results_TMB_LNM), 'signatures')[match(names_sigs_short,
                                                                                                                                 names(python_like_select_name(results_TMB_LNM, 'signatures')))]], give_summary_per_sample)
                         
)
df_roo_sigs = (melt(df_roo_sigs, measure.vars = c('convergence_M', 'convergence_DM', 'convergence_LNM')))
df_roo_sigs$Model = gsub("convergence_", "", df_roo_sigs$variable)
df_roo_sigs[which(df_roo_sigs$num_sigs == 67),]
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
# Full RE intercept
df_roo_sigs_fullRE = data.frame(name=sapply(names_sigs,
                                            function(i) strsplit(i, '_')[[1]][1]),
                                num_samples=as.numeric(sapply(roo_files[sort(names(roo_files[grepl('signatures', names(roo_files))]))],
                                                              function(i) try(nrow(i@count_matrices_all[[1]])))),
                                num_sigs=sapply(roo_files[sort(names(roo_files[grepl('signatures', names(roo_files))]))], get_sigs),
                                convergence_M = sapply(results_TMB_fullRE_M[python_like_select(names(results_TMB_fullRE_M), 'signatures')[match(names_sigs_short,
                                                                                                                                                names(python_like_select_name(results_TMB_fullRE_M, 'signatures')))]], give_summary_per_sample),
                                convergence_DM = sapply(results_TMB_fullRE_DM[python_like_select(names(results_TMB_fullRE_DM), 'signatures')[match(names_sigs_short,
                                                                                                                                                   names(python_like_select_name(results_TMB_fullRE_DM, 'signatures')))]], give_summary_per_sample)
)
df_roo_sigs_fullRE = (melt(df_roo_sigs_fullRE, measure.vars = c('convergence_M', 'convergence_DM')))
df_roo_sigs_fullRE$Model = gsub("convergence_", "", df_roo_sigs_fullRE$variable)
#----------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------#
ggplot()+
  geom_line(data = x, aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  geom_point(data = x, aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  # geom_point(data=df_roo_sigs, aes(y=num_samples, x=1, labels=name))+
  geom_segment(data = df_roo_sigs, aes(y = num_samples, x = 0, xend=num_sigs, yend=num_samples), alpha=0.2, linetype='dashed')+
  geom_segment(data = df_roo_sigs, aes(y = 0, x = num_sigs, yend=num_samples, xend=num_sigs), alpha=0.2, linetype='dashed')+
  geom_label_repel(data = df_roo_sigs, aes(y = num_samples, x = num_sigs, label=name), alpha=0.9, size=4)+
  # geom_label_repel(data = df_roo_sigs, aes(y = 0, x = num_sigs, label=name, direction="y"), alpha=0.2)+
  xlim(c(0, 20))+
  #, ylim=c(0.1, max(df_roo_sigs$num_samples)))+
  scale_y_continuous(trans='log2')+
  labs(x='Number of signatures', y=latex2exp::TeX('Number of samples (log_2)'))+
  ggtitle('Number of signatures and samples in each cancer type, and number of parameters to infer in each model')
# ggsave("../results/assessing_models/num_samples_signatures_parameters.pdf", width = 9, height = 9)

ggplot()+
  geom_line(data = x, aes(x=Var1, y=value, shape=(Effects)))+
  geom_point(data = x, aes(x=Var1, y=value, shape=(Effects)))+facet_wrap(.~factor(Model, levels=c('M', 'DM', 'LNM')))+
  # geom_point(data=df_roo_sigs, aes(y=num_samples, x=1, labels=name))+
  geom_segment(data = df_roo_sigs, aes(y = num_samples, x = 0, xend=num_sigs, yend=num_samples), alpha=0.2, linetype='dashed')+
  geom_segment(data = df_roo_sigs, aes(y = 0, x = num_sigs, yend=num_samples, xend=num_sigs), alpha=0.2, linetype='dashed')+
  geom_label_repel(data = df_roo_sigs, aes(y = num_samples, x = num_sigs, label=name, col=value), alpha=0.9, size=3)+
  # geom_label_repel(data = df_roo_sigs, aes(y = 0, x = num_sigs, label=name, direction="y"), alpha=0.2)+
  xlim(c(0, 20))+
  #, ylim=c(0.1, max(df_roo_sigs$num_samples)))+
  scale_y_continuous(trans='log2')+
  labs(x='Number of signatures', y=latex2exp::TeX('Number of samples (log_2)'))+
  ggtitle('Number of signatures and samples in each cancer type, and number of parameters to infer in each model')
ggsave("../results/assessing_models/num_samples_signatures_parameters_convergence.pdf", width = 16, height = 9)

ggplot()+
  geom_line(data = x, aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  geom_point(data = x, aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  # geom_point(data=df_roo_sigs, aes(y=num_samples, x=1, labels=name))+
  geom_segment(data = df_roo_sigs, aes(y = num_samples, x = 0, xend=num_sigs, yend=num_samples), alpha=0.2, linetype='dashed')+
  geom_segment(data = df_roo_sigs, aes(y = 0, x = num_sigs, yend=num_samples, xend=num_sigs), alpha=0.2, linetype='dashed')+
  geom_label_repel(data = df_roo_sigs, aes(y = num_samples, x = num_sigs, label=name), alpha=0.9, size=4)+
  # geom_label_repel(data = df_roo_sigs, aes(y = 0, x = num_sigs, label=name, direction="y"), alpha=0.2)+
  #, ylim=c(0.1, max(df_roo_sigs$num_samples)))+
  scale_y_continuous(trans='log2')+
  labs(x='Number of signatures', y=latex2exp::TeX('Number of samples (log_2)'))+
  ggtitle('Number of signatures and samples in each cancer type, and number of parameters to infer in each model')

x$Effects[is.na(x$Effects)] = 'Single'
df_roo_sigs_nondup = df_roo_sigs[!duplicated(df_roo_sigs$name),]
df_roo_sigs$Enough_samples = apply(df_roo_sigs, 1, function(ct) DM_ME(as.numeric(ct['num_sigs'])) < as.numeric(ct['num_samples']) )
ggplot()+
  geom_line(data = x[!grepl("LNM", x$Var2),], aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  geom_point(data = x[!grepl("LNM", x$Var2),], aes(x=Var1, y=value, col=(Model), shape=(Effects)))+
  geom_segment(data = df_roo_sigs[!is.na(df_roo_sigs$num_sigs),], aes(y = num_samples, x = 0, xend=num_sigs, yend=num_samples), alpha=0.2, linetype='dashed')+
  geom_segment(data = df_roo_sigs[!is.na(df_roo_sigs$num_sigs),], aes(y = 0, x = num_sigs, yend=num_samples, xend=num_sigs), alpha=0.2, linetype='dashed')+
  geom_label_repel(data = df_roo_sigs_nondup[!is.na(df_roo_sigs_nondup$num_sigs),],
                   aes(y = num_samples, x = num_sigs, label=name),
                   # aes(y = runif(sum(!is.na(df_roo_sigs_nondup$num_sigs)), min = -40, max = 0), x = num_sigs, label=name),
                   alpha=0.9, size=3)+
  xlim(c(0, 20))+
  facet_wrap(.~Enough_samples)+
  # scale_y_continuous(trans='log2')+
  labs(x='Number of signatures', y=latex2exp::TeX('Number of samples (log_2)'))+
  ggtitle('Number of signatures and samples in each cancer type, and number of parameters to infer in each model')
ggsave("../results/results_TMB/pcawg/assessing_models/num_samples_signatures_parameters_split.pdf", width = 16, height = 9)

#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
## Full RE
ggplot()+
  geom_line(data = x[x$Model %in% c('M', 'DM'),], aes(x=Var1, y=value, shape=(Effects)))+
  geom_point(data = x[x$Model %in% c('M', 'DM'),], aes(x=Var1, y=value, shape=(Effects)))+facet_wrap(.~factor(Model, levels=c('M', 'DM')))+
  geom_segment(data = df_roo_sigs_fullRE, aes(y = num_samples, x = 0, xend=num_sigs, yend=num_samples), alpha=0.2, linetype='dashed')+
  geom_segment(data = df_roo_sigs_fullRE, aes(y = 0, x = num_sigs, yend=num_samples, xend=num_sigs), alpha=0.2, linetype='dashed')+
  geom_label_repel(data = df_roo_sigs_fullRE, aes(y = num_samples, x = num_sigs, label=name, col=value), alpha=0.9, size=3)+
  xlim(c(0, 20))+
  scale_y_continuous(trans='log2')+
  labs(x='Number of signatures', y=latex2exp::TeX('Number of samples (log_2)'))+
  ggtitle('Number of signatures and samples in each cancer type, and number of parameters to infer in each model')
ggsave("../results/assessing_models/num_samples_signatures_fullRE_parameters_convergence.pdf", width = 12, height = 9)

ggplot(df_roo_sigs[df_roo_sigs$Model == 'M',], aes(x=num_samples, y=num_sigs, label=name))+
  geom_point()+
  scale_y_continuous(trans='log2')+
  geom_label_repel()
ggsave("../results/assessing_models/num_samples_signatures_scatterplot.pdf", width = 12, height = 9)
#----------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------#
## Looking at particular CT which should have enough samples but which had non-pd results
df_roo_sigs[grep('Lymph-CLL', df_roo_sigs$name),]
results_TMB_fullRE_DM$`Lymph-CLLsignatures`
#----------------------------------------------------------------------------------------#

subset_roo_files = roo_files[!grepl('nucleotidesubstitution3', names(roo_files))]
df_summary = sapply(subset_roo_files, function(i){
  num_active = ncol(attr(i,"count_matrices_active")[[2]])
  if(is.null(num_active)){
    num_active = ""
  }
  c(ct=unique(attr(i,"cancer_type")),
    number_of_samples=nrow(attr(i,"count_matrices_all")[[1]]),
    number_of_categories_1 = ncol(attr(i,"count_matrices_all")[[2]]),
    number_of_categories_2 = num_active)
})

df_summary = do.call('rbind', df_summary)

df_summary_names = do.call('rbind', lapply(gsub("_ROO.RDS", "", basename(rownames(df_summary))),
                                           function(i)  strsplit(i, '_')[[1]]))
# names(df_summary) = basename(names(subset_roo_files))
library(xtable)
print(xtable::xtable(cbind.data.frame(df_summary_names, df_summary[,-1])), include.rownames = FALSE)



