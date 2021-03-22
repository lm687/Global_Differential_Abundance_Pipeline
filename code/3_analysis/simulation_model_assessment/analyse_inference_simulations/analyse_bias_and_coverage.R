rm(list = ls())
require(TMB)
idx_dataset = 1
model <- 'fullREDM'
# model <- 'diagREDM'

fles = list.files("../../../../data/assessing_models_simulation/datasets/", full.names = T)
fles <- fles[grep("multiple_GenerationCnorm_*", fles)]
fles <- fles[grep(paste0("betaintercept", idx_dataset, "_betaslope", idx_dataset, "*"), fles)]
x <- lapply(fles, readRDS)
lst <- list.files("../../../../data/assessing_models_simulation/inference_results/TMB/", full.names = T)
lst <- lst[grepl('multiple_GenerationCnorm_', lst)]
lst <- lst[grepl(model, lst)]
lst <- lst[grep(paste0("betaintercept", idx_dataset, "_betaslope", idx_dataset, "*"), lst)]

all_pd <- lapply(lst, function(i){x <- readRDS(i); try(x$pdHess)})
all_pd[sapply(all_pd, typeof) == 'character'] = FALSE
all_pd_list <- as.vector(unlist(all_pd))

table(all_pd_list)
lst = lst[all_pd_list]
runs <- lapply(lst, readRDS)


summaries = lapply(runs, function(i){
  summary <- summary.sdreport(i)
  summary <- summary[!grepl("u_large", rownames(summary)),]
  # summary[,rownames(summary) == "logs_sd_RE"] = exp(summary[,rownames(summary) == "logs_sd_RE"])
  # summary[,rownames(summary) == "log_lambda"] = exp(summary[,rownames(summary) == "log_lambda"])
  # rownames(summary) = gsub("logs_", "", rownames(summary))
  # rownames(summary) = gsub("log_", "", rownames(summary))
  summary
})
  
summaries_melt = data.frame(melt(sapply(summaries, function(j) j[,1])), stringsAsFactors = F)
summaries_melt$Var1 = as.character(summaries_melt$Var1)
summaries_melt$idx_param = rep(1:sum(summaries_melt$Var2 == 1), length(lst))
summaries_melt[summaries_melt$Var1 == "log_lambda", "value"] = exp(summaries_melt[summaries_melt$Var1 == "log_lambda", "value"])
summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"] = exp(summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "value"])

summaries_melt[summaries_melt$Var1 == "log_lambda", "Var1"] = "lambda"
summaries_melt[summaries_melt$Var1 == "logs_sd_RE", "Var1"] = "sd_RE"


if(model == "fullREDM"){
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                rep(0,((x[[1]]$d-1)**2-(x[[1]]$d-1))/2), ## covariances RE
                rep(x[[1]]$sd_RE, x[[1]]$d-1), ##sd RE ## this is particular to GenerationCnorm
                x[[1]]$lambda)
}else if(model == "diagREDM"){
  true_vals = c(as.vector(x[[1]]$beta), ## betas
                rep(x[[1]]$sd_RE, x[[1]]$d-1), ##sd RE ## this is particular to GenerationCnorm
                x[[1]]$lambda)
}
summaries_melt$subtract = summaries_melt$value - rep(true_vals, length(lst))

ggplot(summaries_melt, aes(x=idx_param, y=subtract, group=idx_param))+
  geom_abline(slope = 0, intercept = 0, lty='dashed', col='blue')+
  geom_boxplot()+facet_wrap(~Var1, scales = "free", nrow=1)+labs(x="Parameter", y="Bias")
ggsave(paste0("~/Desktop/boxplots/setsim_",model, '_', idx_dataset, "_bias.pdf"), width = 10, height = 3.5)
# summaries_mat = matrix(summaries_melt$value, nrow=sum(summaries_melt$Var2 == 1))


(x[[1]]$beta)[1,1]
sapply(summaries, function(j) j[1,1])

(x[[1]]$beta)[2,1]
sapply(summaries, function(j) j[2,1])

boxplot(sapply(summaries, function(j) j[2,1]))
abline(h = (x[[1]]$beta)[2,1])

x[[1]]$lambda
exp(sapply(summaries, function(j) j[rownames(j) == "log_lambda",1]))

## compute coverage
## for each estimate, check if the true value is in the 95% confidence interval computed from the standard errors
summaries[[1]]
true_vals
confints <- sapply(summaries, function(it_run){
sapply(1:length(true_vals), function(i){
  confintint = c(it_run[i,1]-1.96*it_run[i,2],it_run[i,1]+1.96*it_run[i,2])
  (true_vals[i] > confintint[1]) & (true_vals[i] < confintint[2])
  } )
})
rownames(confints) = rownames(summaries[[1]])

ggplot(cbind.data.frame(parameter=make.names(rownames(confints), unique = T), CI=apply(confints, 1, mean)),
       aes(x=parameter, y=CI, group=1))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("~/Desktop/boxplots/setsim_",model, '_', idx_dataset, "_coverage.pdf"), width = 10, height = 3.5)


