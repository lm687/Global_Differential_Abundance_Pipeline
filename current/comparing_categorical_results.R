#------------------------------------------------------------------------------------------------------------------#
rm(list = ls())

setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/")

source("../1_create_ROO/roo_functions.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")

library(TMB)
library(dplyr)

TMB::compile("mm_multinomial/fullRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial"))

TMB::compile("mm_multinomial/fullRE_ME_multinomial_categorical.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial_categorical"))
#------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------#
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                                   strsplit, split = "_")))
colnames(samples_files) = c('CT', 'type')
# table(samples_files[,1], samples_files[,2])
# ct = "Bladder-TCC" #samples_files[1,1]
# typedata =nucleotidesubstitution3  #"signatures" #samples_files[1,2]

samples_files2 = samples_files %>% filter(type != "nucleotidesubstitution3")
rownames(samples_files2) = rownames(samples_files)[samples_files$type != "nucleotidesubstitution3"]
#------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------#
## get some cancer type and compare multinomial and categorical results
results_m = readRDS("../../data/pcawg_robjects_cache/tmb_results/fullRE_M_Bladder-TCC_signatures.RDS")
results_m

image(results_m$cov.fixed)
image(results_m$cov.fixed[c(T,F),c(T,F)])
image(results_m$cov.fixed[c(F,T),c(F,T)])

idx = 4

i = samples_files2[idx,]
i
# results_Mcat = (wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_Mcat"))
# results_Mcat
# 
# plot(results_m$par.fixed,
#      results_Mcat$par.fixed)
#------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------#
folder_robjs = "../../data/pcawg_robjects_cache/tmb_results/"
results_TMB_fullRE_Mcat = lapply( python_like_select(list.files(folder_robjs), "^fullRE_Mcat_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_Mcat) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_Mcat_"), clean_name_fullRE)

results_TMB_fullRE_DMcat = lapply( python_like_select(list.files(folder_robjs), "^fullRE_DMcat_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DMcat) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_DMcat_"), clean_name_fullRE)

results_TMB_fullRE_M = lapply( python_like_select(list.files(folder_robjs), "^fullRE_M_"), function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_M) = sapply(python_like_select(list.files(folder_robjs), "^fullRE_M_"), clean_name_fullRE)

full_RE_DM = python_like_select(list.files(folder_robjs), "^fullRE_DM_"); full_RE_DM = full_RE_DM[-grep("_altpar_", full_RE_DM)]
results_TMB_fullRE_DM = lapply( full_RE_DM, function(i) readRDS(paste0(folder_robjs, i)))
names(results_TMB_fullRE_DM) = sapply(full_RE_DM, clean_name_fullRE)
#------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------#
# results_TMB_fullRE_M = results_m

results_TMB_fullRE_Mcat
results_TMB_fullRE_M

all(names(results_TMB_fullRE_Mcat) == names(results_TMB_fullRE_M))

## in general there is agreement on whwther the model converged
table(good_M=unlist(sapply(results_TMB_fullRE_M, `[`, "pdHess")),
      good_Mcat=unlist(sapply(results_TMB_fullRE_Mcat, `[`, "pdHess")))

## select cases where both converged, only
which_good_conv = as.logical(unlist(sapply(results_TMB_fullRE_M, `[`, "pdHess"))) & as.logical(unlist(sapply(results_TMB_fullRE_Mcat, `[`, "pdHess")))

## there is good agreement for the estimates
plot(log(unlist(sapply(results_TMB_fullRE_Mcat[which_good_conv], `[`, "par.fixed"))),
       log(unlist(sapply(results_TMB_fullRE_M[which_good_conv], `[`, "par.fixed"))))

### same for DM

length(results_TMB_fullRE_DMcat)
length(results_TMB_fullRE_DM)
which_good_conv_M = as.logical(unlist(sapply(results_TMB_fullRE_M, `[`, "pdHess"))) & as.logical(unlist(sapply(results_TMB_fullRE_Mcat, `[`, "pdHess")))
which_good_conv_DM = as.logical(unlist(sapply(results_TMB_fullRE_DM, `[`, "pdHess"))) & as.logical(unlist(sapply(results_TMB_fullRE_DMcat, `[`, "pdHess")))

results_TMB_fullRE_DM_subset = results_TMB_fullRE_DM[match(names(results_TMB_fullRE_DMcat),names(results_TMB_fullRE_DM))]
results_TMB_fullRE_M_subset = results_TMB_fullRE_M[match(names(results_TMB_fullRE_Mcat),names(results_TMB_fullRE_M))]

# which_good_conv_M = ( ((sapply(sapply(results_TMB_fullRE_M_subset, `[`, "par.fixed"), typeof) ) == "double") &
#                          ((sapply(sapply(results_TMB_fullRE_Mcat, `[`, "par.fixed"), typeof) ) == "double"))
# which_good_conv_DM = ( ((sapply(sapply(results_TMB_fullRE_DM_subset, `[`, "par.fixed"), typeof) ) == "double") &
#                       ((sapply(sapply(results_TMB_fullRE_DMcat, `[`, "par.fixed"), typeof) ) == "double"))

plot(log(unlist(sapply(results_TMB_fullRE_DM_subset[which_good_conv_DM], `[`, "par.fixed"))),
     log(unlist(sapply(results_TMB_fullRE_DMcat[which_good_conv_DM], `[`, "par.fixed"))))
abline(coef=c(0,1), lty='dashed')

require(ggrepel)
df_m = cbind.data.frame(M=(unlist(sapply(results_TMB_fullRE_M_subset[which_good_conv_M], `[`, "par.fixed"))),
                         Mcat=(unlist(sapply(results_TMB_fullRE_Mcat[which_good_conv_M], `[`, "par.fixed"))),
                         parameter_name = unlist(sapply(results_TMB_fullRE_M_subset[which_good_conv_M], function(i) names(i$par.fixed))),
                         dataset=rep(names(results_TMB_fullRE_M_subset[which_good_conv_M]), sapply(results_TMB_fullRE_M_subset[which_good_conv_M], function(i) length(i$par.fixed))))
df_m = cbind.data.frame(M=log(unlist(sapply(results_TMB_fullRE_M_subset[which_good_conv_M], `[`, "par.fixed"))),
                        Mcat=log(unlist(sapply(results_TMB_fullRE_Mcat[which_good_conv_M], `[`, "par.fixed"))),
                        parameter_name = unlist(sapply(results_TMB_fullRE_M_subset[which_good_conv_M], function(i) names(i$par.fixed))),
                        dataset=rep(names(results_TMB_fullRE_M_subset[which_good_conv_M]), sapply(results_TMB_fullRE_M_subset[which_good_conv_M], function(i) length(i$par.fixed))))
df_dm = cbind.data.frame(DM=(unlist(sapply(results_TMB_fullRE_DM_subset[which_good_conv_DM], `[`, "par.fixed"))),
                 DMcat=(unlist(sapply(results_TMB_fullRE_DMcat[which_good_conv_DM], `[`, "par.fixed"))),
                 parameter_name = unlist(sapply(results_TMB_fullRE_DM_subset[which_good_conv_DM], function(i) names(i$par.fixed))),
                 dataset=rep(names(results_TMB_fullRE_DM_subset[which_good_conv_DM]), sapply(results_TMB_fullRE_DM_subset[which_good_conv_DM], function(i) length(i$par.fixed))))

#' If including non-converged: In the comparison between M and Mcat, there are clearly two groups of CT: those for which there is a striking
#' agreement in all parameters, and those for which there is a substantially worse agreement in all parameters.
#' It seems as though this bad agreement is in the signaturesm, and not nucleotide subs., and in cancer types in which there is a very 
#' high number of active signatures.
#' If excluding non-converged: perfect agreement
ggplot(droplevels(df_m),
       aes(x=M, y=Mcat, col=parameter_name, label=dataset))+geom_point()+#guides(col=FALSE)+
  geom_abline(slope = 1, intercept = 0)+facet_wrap(.~dataset, scales = "free")
ggsave("~/Desktop/M_vs_catM_scatter.png", width = 10, height = 10)

#' If including non-converged: There isn't any agreement between the dirichlet and the categorical dirichlet. The values of beta are always much lower in the 
#' catgeorical case. Covariances and variances are also not particularly in agreement, although closer
#' If excluding non-converged: pretty good convergence
ggplot(df_dm,
       aes(x=DM, y=DMcat, col=parameter_name, label=dataset))+geom_point()+#guides(col=FALSE)+
  geom_abline(slope = 1, intercept = 0)+facet_wrap(.~dataset)
ggsave("~/Desktop/DM_vs_catDM_scatter.png", width = 10, height = 10)
#------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------#
df_joint = cbind.data.frame(DM=df_dm[df_dm$parameter_name != "log_lambda",],
      M=df_m[df_m$dataset %in% intersect(unique(df_dm$dataset), unique(df_m$dataset)),])
all(as.character(df_joint$DM.dataset) == as.character(df_joint$M.dataset))

df_joint$M.

##' In the M vs DM comparison, beta is generally in very good agreement (with a very few exceptions)
##' Of course, few samples managed to converge in the DM case
##' The parameters for the random effects can vary substantially between M and DM
ggplot(df_joint,
       aes(x=M.M, y=DM.DM, col=M.parameter_name, label=M.dataset))+geom_point()+#guides(col=FALSE)+
geom_abline(slope = 1, intercept = 0)+facet_wrap(.~M.dataset)
#------------------------------------------------------------------------------------------------------------------#





