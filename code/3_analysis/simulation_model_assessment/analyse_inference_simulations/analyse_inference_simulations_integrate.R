rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

runs_M = readRDS("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/GenerationCnorm_fullREM.RDS")
runs_DM = readRDS("../../../../data/assessing_models_simulation/inference_results/TMB/summaries/GenerationCnorm_fullREDM.RDS")

joint_df = cbind.data.frame(M=runs_M, DM=runs_DM[match(rownames(runs_M), rownames(runs_DM)),])

pairs(cbind(beta_true=joint_df$M.beta_true, M_beta_est=joint_df$M.beta_est, DM_beta_est=joint_df$DM.beta_est))

pdf("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/GenerationCnorm/summaries/M_DM_comparison.pdf")
do.call('grid.arrange', list(ggplot(joint_df, aes(x=M.beta_true, M.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1),
                  ggplot(joint_df, aes(x=DM.beta_true, DM.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)))
dev.off()

pdf("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/GenerationCnorm/summaries/M_DM_comparison_only_common.pdf")
do.call('grid.arrange', list(ggplot(joint_df[!is.na(joint_df$M.beta_est) & !is.na(joint_df$DM.beta_est),], aes(x=M.beta_true, M.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1),
                             ggplot(joint_df[!is.na(joint_df$M.beta_est) & !is.na(joint_df$DM.beta_est),], aes(x=DM.beta_true, DM.beta_est))+geom_point()+geom_abline(intercept = 0, slope = 1)))
dev.off()



runs_ttest_irl = lapply(datasets_files, function(i)  try(wrapper_run_ttest_ilr(i)))
runs_ttest_props = lapply(datasets_files, function(i)  try(wrapper_run_ttest_props(i)))
pvals_ttest_ilr = as.numeric(unlist(runs_ttest_irl))
pvals_ttest_ilr_adj = pvals_ttest_ilr

res_M = readRDS("../data/assessing_models_simulation/inference_results/TMB/summaries/GenerationCnorm_fullREM.RDS")

pvals_m = res_M[sapply(unique(res_M$idx),
                       function(i) which(res_M$idx == i)[1]),'pvals_adj']
res_all = rbind(summarise_DA_detection(true = DA_bool, predicted = pvals_m < 0.05),
                summarise_DA_detection(true = DA_bool, predicted = pvals_adj <= 0.05),
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
ggsave("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/GenerationCnorm/summaries/M_betaslopes_confint.pdf",
       height = 14, width = 8)

ggplot(joint_df,
       aes(x=DM.idx_within_dataset, y=(DM.beta_est), col=DM.pvals_adj<0.05))+
  geom_point(aes(x=M.idx_within_dataset, y=M.beta_true), shape=4)+
  geom_abline(slope = 0, intercept = 0, alpha=0.2)+geom_point(aes(shape=DM.DA_bool))+
  geom_errorbar(aes(ymin=DM.beta_est-1.96*DM.beta_stderr, ymax=DM.beta_est+1.96*DM.beta_stderr), width=.2,
                position=position_dodge(.9))+
  facet_wrap(.~interaction(DM.idx, DM.beta_gamma_shape, DM.DA_bool), scales='free_x', nrow=length(unique(joint_df$DM.beta_gamma_shape)))+
  theme_bw()+theme(legend.position = "bottom")
ggsave("../../../../results/results_TMB/simulated_datasets/mixed_effects_models/GenerationCnorm/summaries/DM_betaslopes_confint.pdf",
       height = 14, width = 8)

## What hasn't run?
cat(paste0(sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) substr(i, start = 1, stop = nchar(i)-1)) %>% unique, sep='" "', collapse=''))

sapply(rownames(joint_df[is.na(joint_df$DM.beta_est),]), function(i) gsub("_dataset", "", substr(i, start = 1, stop = nchar(i)-1))) %>% unique

all(joint_df$M.DA_bool == joint_df$DM.DA_bool)
table(truth=joint_df$M.DA_bool, M=joint_df$M.pvals_adj <= 0.05)
table(truth=joint_df$DM.DA_bool, DM=joint_df$DM.pvals_adj <= 0.05)

