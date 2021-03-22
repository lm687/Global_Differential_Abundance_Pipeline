rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../../2_inference_TMB/helper_TMB.R")

library(gridExtra)

generation = "GenerationCnorm"
generation = "generationFnorm"
generation = "generationGnorm"
generation = "generationMGnorm"

runs_fullREM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", generation, "_fullREM.RDS"))
runs_fullREDM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", generation, "_fullREDM.RDS"))

## Problem with convergence is acute in DM
table(is.na(runs_fullREM$beta_est))
table(is.na(runs_fullREDM$beta_est))

joint_df = cbind.data.frame(M=runs_fullREM,
                            DM=runs_fullREDM[match(rownames(runs_fullREM),
                                                                   rownames(runs_fullREDM)),])

pairs(cbind(beta_true=joint_df$M.beta_true, M_beta_est=joint_df$M.beta_est, DM_beta_est=joint_df$DM.beta_est))

system(paste0("mkdir -p ../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/"))

pdf(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M_DM_comparison.pdf"))
do.call('grid.arrange', list(ggplot(joint_df, aes(x=M.beta_true, M.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1),
                  ggplot(joint_df, aes(x=DM.beta_true, DM.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)))
dev.off()

pdf(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M_DM_comparison_only_common.pdf"))
do.call('grid.arrange', list(ggplot(joint_df[!is.na(joint_df$M.beta_est) & !is.na(joint_df$DM.beta_est),], aes(x=M.beta_true, M.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1),
                             ggplot(joint_df[!is.na(joint_df$M.beta_est) & !is.na(joint_df$DM.beta_est),], aes(x=DM.beta_true, DM.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)))
dev.off()

datasets_files = list.files("../../../../data/assessing_models_simulation/datasets/", full.names = TRUE)
datasets_files = datasets_files[grep(pattern = generation, datasets_files)]

# match
datasets_files = datasets_files[match(unique(sapply(rownames(joint_df), function(i) strsplit(i, "_dataset")[[1]][1])),
                                      gsub("_dataset.RDS", "", basename(datasets_files)))]

datasets = lapply(datasets_files, readRDS)
names(datasets) = basename(datasets_files)
DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )

runs_ttest_irl = lapply(datasets_files, function(i)  try(wrapper_run_ttest_ilr(i)))
runs_ttest_props = lapply(datasets_files, function(i)  try(wrapper_run_ttest_props(i)))
pvals_ttest_ilr = as.numeric(unlist(runs_ttest_irl))
pvals_ttest_ilr_adj = pvals_ttest_ilr

res_M = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", generation, "_fullREM.RDS"))
res_DM = readRDS(paste0("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/", generation, "_fullREDM.RDS"))

pvals_m = res_M[sapply(unique(res_M$idx),
                       function(i) which(res_M$idx == i)[1]),'pvals_adj']
pvals_dm = res_DM[sapply(unique(res_DM$idx),
                       function(i) which(res_DM$idx == i)[1]),'pvals_adj']
res_all = rbind(summarise_DA_detection(true = DA_bool, predicted = pvals_m < 0.05),
                summarise_DA_detection(true = DA_bool, predicted = pvals_dm <= 0.05),
                # summarise_DA_detection(true = DA_bool, predicted = runs_ttest_props <= 0.05),
                summarise_DA_detection(true = DA_bool, predicted = pvals_ttest_ilr_adj <= 0.05))
rownames(res_all) = c('Multinomial', 'Dirichlet-Multinomial', 'ILR')
res_all
xtable::xtable(res_all)

table(DA_bool, M_est=pvals_m <= 0.05)
table(DA_bool, DM_est=pvals_adj <= 0.05)
table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05)

xtable::xtable(table(DA_bool, M_est=pvals_m <= 0.05))
xtable::xtable(table(DA_bool, DM_est=pvals_adj <= 0.05))
xtable::xtable(table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05))

require(reshape2)
head(melt(list(table(DA_bool, M_est=pvals_m <= 0.05),
               table(DA_bool, DM_est=pvals_adj <= 0.05),
               table(DA_bool, ILR_est=pvals_ttest_ilr_adj <= 0.05))
))


joint_df$M.beta_true


ggplot(joint_df,
       aes(x=M.idx_within_dataset, y=(M.beta_est), col=M.pvals_adj<0.05))+
  geom_point(aes(x=M.idx_within_dataset, y=M.beta_true), shape=4)+
  geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point(aes(shape=M.DA_bool))+
  geom_errorbar(aes(ymin=M.beta_est-1.96*M.beta_stderr, ymax=M.beta_est+1.96*M.beta_stderr), width=.2,
                position=position_dodge(.9))+
  facet_wrap(.~interaction(M.idx, M.beta_gamma_shape, M.DA_bool), scales='free_x', nrow=length(unique(joint_df$M.beta_gamma_shape)))+
  theme_bw()+theme(legend.position = "bottom")
  # scale_colour_viridis_d(option = "plasma")
  # scale_colour_manual(values = c("red","#2e8b57", "red", "#2e8b57")) #"#2e8b57"))
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M_betaslopes_confint.pdf"),
       height = 14, width = 8)

ggplot(joint_df,
       aes(x=DM.idx_within_dataset, y=(DM.beta_est), col=DM.pvals_adj<0.05))+
  geom_point(aes(x=M.idx_within_dataset, y=M.beta_true), shape=4)+
  geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point(aes(shape=DM.DA_bool))+
  geom_errorbar(aes(ymin=DM.beta_est-1.96*DM.beta_stderr, ymax=DM.beta_est+1.96*DM.beta_stderr), width=.2,
                position=position_dodge(.9))+
  facet_wrap(.~interaction(DM.idx, DM.beta_gamma_shape, DM.DA_bool), scales='free_x', nrow=length(unique(joint_df$DM.beta_gamma_shape)))+
  theme_bw()+theme(legend.position = "bottom")
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/DM_betaslopes_confint.pdf"),
       height = 14, width = 8)

all(joint_df$M.DA_bool == joint_df$DM.DA_bool)
table(truth=joint_df$M.DA_bool, M=joint_df$M.pvals_adj <= 0.05)
table(truth=joint_df$DM.DA_bool, DM=joint_df$DM.pvals_adj <= 0.05)

## ROC curves

colnames(joint_df)

all(joint_df$M.d == joint_df$DM.d); all(joint_df$M.n == joint_df$DM.n); all(joint_df$M.beta_gamma_shape == joint_df$DM.beta_gamma_shape);  all(joint_df$M.DA_bool == joint_df$M.DA_bool); 
# "M.d", "M.n", "M.beta_gamma_shape", M.pvals_adj", "M.DA_bool", "DM.pvals_adj" 
## select datasets, not betas

joint_df$M_type = paste0(joint_df$DM.DA_bool, joint_df$M.pvals_adj < 0.05)
joint_df$DM_type = paste0(joint_df$DM.DA_bool, joint_df$DM.pvals_adj < 0.05)
joint_df_grouping_by_n = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M.DA_bool", "DM.pvals_adj", "M_type", "DM_type" )) %>% 
  group_by(M.beta_gamma_shape, M.n) %>%
  mutate(sensitivity_DM=sum(DM_type == 'TRUETRUE')/sum(DM_type %in% c('TRUETRUE', 'TRUEFALSE')),
         sensitivity_M=sum(M_type == 'TRUETRUE')/sum(M_type %in% c('TRUETRUE', 'TRUEFALSE')),
         specificity_DM=sum(DM_type == 'FALSEFALSE')/sum(DM_type %in% c('FALSETRUE', 'FALSEFALSE')),
         specificity_M=sum(M_type == 'FALSEFALSE')/sum(M_type %in% c('FALSETRUE', 'FALSEFALSE')))

joint_df_grouping_by_d = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M.DA_bool", "DM.pvals_adj", "M_type", "DM_type" )) %>% 
  group_by(M.beta_gamma_shape, M.d) %>%
  mutate(sensitivity_DM=sum(DM_type == 'TRUETRUE')/sum(DM_type %in% c('TRUETRUE', 'TRUEFALSE')),
         sensitivity_M=sum(M_type == 'TRUETRUE')/sum(M_type %in% c('TRUETRUE', 'TRUEFALSE')),
         specificity_DM=sum(DM_type == 'FALSEFALSE')/sum(DM_type %in% c('FALSETRUE', 'FALSEFALSE')),
         specificity_M=sum(M_type == 'FALSEFALSE')/sum(M_type %in% c('FALSETRUE', 'FALSEFALSE')))


ggplot(droplevels(joint_df_grouping_by_n), aes(x=sensitivity_M, y=1-specificity_M))+
  geom_point()
ggplot(droplevels(joint_df_grouping_by_d), aes(x=sensitivity_M, y=1-specificity_M))+
  geom_point()

ggplot(droplevels(joint_df_grouping_by_n), aes(x=sensitivity_DM, y=1-specificity_DM, col=interaction(M.d, M.n)))+
  geom_point()

ggplot(droplevels(joint_df_grouping_by_n), aes(x=M.n, y=sensitivity_DM, group=M.beta_gamma_shape))+
  geom_point()+geom_line()+facet_wrap(.~M.beta_gamma_shape)
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M.beta_gamma_shape_sensitivity.pdf"))

head(melt(joint_df_grouping_by_n, id.vars = c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M_type", "M.DA_bool", "sensitivity_M", "specificity_DM")))
plot(joint_df_grouping_by_n$sensitivity_M, 1-joint_df_grouping_by_n$specificity_M)

head(melt(joint_df_grouping_by_n[,c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M_type", "M.DA_bool", "sensitivity_M", "specificity_DM")],
          id.vars = c("M.d", "M.n", "M.beta_gamma_shape")))

modify_colnames = function(i){
  colnames(i) = gsub("DM.", "", colnames(i))
  colnames(i) = gsub("M.", "", colnames(i))
  colnames(i) = gsub("_M", "", colnames(i))
  colnames(i) = gsub("_DM", "", colnames(i))
  i
}

joint_df_grouping_by_n_alt = rbind(modify_colnames(cbind(joint_df_grouping_by_n[,c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M_type", "M.DA_bool", "sensitivity_M", "specificity_M")], model="M")),
      modify_colnames(cbind(joint_df_grouping_by_n[,c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "DM_type", "M.DA_bool", "sensitivity_DM", "specificity_DM")], model="DM")))

joint_df_grouping_by_d_alt = rbind(modify_colnames(cbind(joint_df_grouping_by_d[,c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "M_type", "M.DA_bool", "sensitivity_M", "specificity_M")], model="M")),
                                   modify_colnames(cbind(joint_df_grouping_by_d[,c("M.d", "M.n", "M.beta_gamma_shape", "M.pvals_adj", "DM_type", "M.DA_bool", "sensitivity_DM", "specificity_DM")], model="DM")))

ggplot(droplevels(joint_df_grouping_by_n_alt), aes(x=d, y=sensitivity, group=interaction(n, model), col=n, shape=model))+
  geom_point()+geom_line()+facet_wrap(.~beta_gamma_shape)

ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=d, y=sensitivity, group=interaction(n, model), col=n, shape=model))+
  geom_point()+geom_line()+facet_wrap(.~beta_gamma_shape)
ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=sensitivity, group=interaction(d, model), col=n, shape=model))+
  geom_point()+geom_line()

ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=sensitivity, group=interaction(beta_gamma_shape, model), col=model))+
  geom_point()+geom_line()+geom_boxplot()+scale_x_continuous(trans = "log2")

ggplot(droplevels(joint_df_grouping_by_d_alt), aes(x=beta_gamma_shape, y=specificity, group=interaction(beta_gamma_shape, model), col=model))+
  geom_point()+geom_line()+geom_boxplot()+scale_x_continuous(trans = "log2")

## What hasn't run?
cat(paste0(sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) substr(i, start = 1, stop = nchar(i)-1)) %>% unique, sep='" "', collapse=''))

sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) gsub("_dataset", "", substr(i, start = 1, stop = nchar(i)-1))) %>% unique

ggplot(joint_df, aes(x=M.d, y=as.numeric(M.converged), group=interaction(M.n, M.d)))+
  geom_jitter(height = 0.1, alpha=0.2, col='#07367d')+geom_violin()+facet_wrap(.~M.n)+
  theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  ggtitle('Success in convergence of Multinomial runs')
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M_d_n_convergence.pdf"))

ggplot(joint_df, aes(x=DM.d, y=as.numeric(DM.converged), group=interaction(DM.n, DM.d)))+
  geom_jitter(height = 0.1, alpha=0.2, col='#07367d')+geom_violin()+facet_wrap(.~DM.n)+
  theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  ggtitle('Success in convergence of Dirichlet-Multinomial runs')
ggsave(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/DM_d_n_convergence.pdf"))

ggplot(joint_df, aes(x=M.d, y=as.numeric(M.converged), group=interaction(M.n, M.d), col=M.beta_gamma_shape))+
  geom_jitter(height = 0.1, alpha=0.8)+geom_violin()+facet_wrap(.~M.n)+
  theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  ggtitle('Success in convergence of Multinomial runs')

ggplot(joint_df, aes(x=M.beta_gamma_shape, y=as.numeric(M.converged), group=M.beta_gamma_shape, col=M.d))+
  geom_violin()+
  geom_point()+
  # geom_jitter(height = 0.1, alpha=0.8)+
  facet_wrap(.~M.n)+
  theme_bw()+labs(x='Number of categories (d)', y='Number of successful (1) or unsuccessful (0) convergences')+
  ggtitle('Success in convergence of Multinomial runs')

## get percentage of successful runs
joint_df_grouping_convergence = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "DM.converged", "M.converged")) %>% 
  group_by(M.n, M.d) %>%
  mutate(convergence_M=sum(M.converged)/(sum(M.converged)+sum(!M.converged)),
         convergence_DM=sum(DM.converged)/length(DM.converged))
joint_df_grouping_convergence_2 = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "DM.converged", "M.converged", "M.beta_gamma_shape")) %>% 
  group_by(M.beta_gamma_shape) %>%
  mutate(convergence_M=sum(M.converged)/(sum(M.converged)+sum(!M.converged)),
         convergence_DM=sum(DM.converged)/length(DM.converged))
joint_df_grouping_convergence_2_n = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "DM.converged", "M.converged", "M.beta_gamma_shape")) %>% 
  group_by(M.beta_gamma_shape, M.n) %>%
  mutate(convergence_M=sum(M.converged)/(sum(M.converged)+sum(!M.converged)),
         convergence_DM=sum(DM.converged)/length(DM.converged))
joint_df_grouping_convergence_2_d = joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>%
  select(c("M.d", "M.n", "DM.converged", "M.converged", "M.beta_gamma_shape")) %>% 
  group_by(M.beta_gamma_shape, M.d) %>%
  mutate(convergence_M=sum(M.converged)/(sum(M.converged)+sum(!M.converged)),
         convergence_DM=sum(DM.converged)/length(DM.converged))

pdf(paste0("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/", generation, "/summaries/M_DM_convergence.pdf"))
grid.arrange(ggplot(joint_df_grouping_convergence, aes(x=M.d, y=convergence_M, col=M.n, group=M.n))+geom_line()+geom_point()+ggtitle('Convergence Multinomial'),
ggplot(joint_df_grouping_convergence, aes(x=M.d, y=convergence_DM, col=M.n, group=M.n))+geom_line()+geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
dev.off()

joint_df[sapply(unique(joint_df$DM.idx), function(i) which(joint_df$DM.idx == i)[1]),] %>% filter(M.d==4, M.n==100) %>%
  select(c("M.d", "M.n", "DM.converged", "M.converged")) 
head(joint_df_grouping_convergence)

grid.arrange(ggplot(joint_df_grouping_convergence_2, aes(x=M.beta_gamma_shape, y=convergence_M, group=M.beta_gamma_shape))+geom_line()+
               geom_point()+ggtitle('Convergence Multinomial'),
             ggplot(joint_df_grouping_convergence_2, aes(x=M.beta_gamma_shape, y=convergence_DM, group=M.beta_gamma_shape))+geom_line()+
               geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
grid.arrange(ggplot(joint_df_grouping_convergence_2_n, aes(x=M.beta_gamma_shape, y=convergence_M, col=M.n, group=interaction(M.n)))+geom_line()+
               geom_point()+ggtitle('Convergence Multinomial'),
             ggplot(joint_df_grouping_convergence_2_n, aes(x=M.beta_gamma_shape, y=convergence_DM, col=M.n,  group=interaction(M.n)))+geom_line()+
               geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))
grid.arrange(ggplot(joint_df_grouping_convergence_2_d, aes(x=M.beta_gamma_shape, y=convergence_M, col=M.d, group=interaction(M.d)))+geom_line()+
               geom_point()+ggtitle('Convergence Multinomial'),
             ggplot(joint_df_grouping_convergence_2_d, aes(x=M.beta_gamma_shape, y=convergence_DM, col=M.d,  group=interaction(M.d)))+geom_line()+
               geom_point()+ggtitle('Convergence Dirichlet-Multinomial'))

createbarplot_object("../../../../data/assessing_models_simulation/datasets//GenerationCnorm_50_100_80_4_0.8_dataset.RDS")
datasets$GenerationCnorm_50_100_80_4_0.8_dataset.RDS$objects_counts
joint_df[grepl('GenerationCnorm_50_100_80_4_0.8_dataset', rownames(joint_df)),]
