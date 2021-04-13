### Troubleshooting inference with TMB
rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/")
library(TMB)
library(ggplot2)
library(dplyr)
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("../../../CDA_in_Cancer/code/functions/meretricious/pretty_plots/prettySignatures.R")

TMB::compile("mm_multinomial/fullRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial"))
TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial"))
TMB::compile("mm_multinomial/diagRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/diagRE_ME_dirichletmultinomial"))
TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial_sparsecov.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial_sparsecov"))

ct = "Kidney-RCC.clearcell" #samples_files[1,1]
typedata =  "signatures" #samples_files[1,2]
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                                   strsplit, split = "_")))
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# results_Mcat = (wrapper_run_TMB(i[1,1], i[1,2], model = "fullRE_Mcat"))
results_M = (wrapper_run_TMB(ct, typedata, model = "diagRE_M"))
results_M

## subset: 1. create subset
.xxx = readRDS(file = paste0("../../data/roo/", ct, '_', typedata, "_ROO.RDS" ))
saveRDS(give_subset_sigs(.xxx, sigs_to_remove = c('SBS3', 'SBS24', 'SBS26')), file = paste0("../../data/roo/", ct, '-subset_', typedata, "_ROO.RDS" ))
saveRDS(give_subset_sigs(.xxx, sigs_to_remove = c('SBS3', 'SBS22', 'SBS24', 'SBS26')), file = paste0("../../data/roo/", ct, '-subset_', typedata, "_ROO.RDS" ))

## subset: 2. infer subset
results_M_subset = (wrapper_run_TMB(paste0(ct, '-subset'), typedata, model = "fullRE_M", sort_columns = T))
results_M_subset
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
## Plot barplot

give_barplot(ct, typedata, simulation = F, legend_on = T)
obj = load_PCAWG(ct, typedata, simulation=F)
View(obj$Y)
dim(obj$Y)
colnames(obj$Y)
obj_mod = obj
# obj_mod$Y = obj_mod$Y[,c(1:9, 12)]
results_M_mod = wrapper_run_TMB(object = obj_mod, model="fullRE_M")
results_M_mod ## this does converge well

## correlation between absolute exposures
results_DMdiag_mod = wrapper_run_TMB(object = obj_mod, model="diagRE_DM")

pairs(log(obj$Y))
plot(obj$Y[,'SBS6'],
     obj$Y[,'SBS12'])

remove_inf = function(i){
  i = i[!is.infinite(i)]
  i[!is.nan(i)]
}

give_LRchanges_barplot = function(obj_it, nrow_facets=NULL, ct_name, typedata){
  mat_LRchanges = melt(log(normalise_rw(obj_it$Y[obj_it$x[,2] == 0,])/normalise_rw(obj_it$Y[obj_it$x[,2] == 1,])))
  mat_LRchanges$col = 'NonInf'
  mat_LRchanges$col[is.infinite(mat_LRchanges$value)] = 'Inf'
  mat_LRchanges$value[is.infinite(mat_LRchanges$value)] = min( remove_inf(mat_LRchanges$value))
  mat_LRchanges$col[is.nan(mat_LRchanges$value)] = 'Inf'
  mat_LRchanges$value[is.nan(mat_LRchanges$value)] = max( remove_inf(mat_LRchanges$value))
  
  a = ggplot(mat_LRchanges,
             aes(x=interaction(Var1, Var2), y=value, fill=col))+geom_bar(stat = 'identity')+
    labs(x='Patients', y='logR change')+
    scale_fill_manual(values=c('Inf'='blue', 'NonInf'='black'))+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.title = element_blank())+
    ggtitle(paste0('Log ratio changes: ', ct_name, ' ', typedata))
  if(is.null(nrow_facets)){
    a = a + facet_wrap(.~Var2, scales = "free_x")
  }else{
    a = a + facet_wrap(.~Var2, scales = "free_x", nrow=nrow_facets)
  }
  return(a)
}

pdf("~/Desktop/logR.pdf")
for(ct_it in sort(unique(samples_files$X1))){
  for(typedata_it in c('nucleotidesubstitution1', 'signatures')){
    obj_it = load_PCAWG(ct_it, typedata_it, simulation=F)

    if(!is.na(obj_it)){
      
      a = give_LRchanges_barplot(obj_it, ct_name=ct_it)
      print(a)
    }
  }
}
dev.off()

obj_breast = load_PCAWG("Breast-AdenoCA", "signatures", simulation=F)
breastadenoca = give_LRchanges_barplot(obj_it = obj_breast,
                                       nrow_facets = 1, ct_name = "Breast-AdenoCA", typedata = "signatures")
breastadenoca
ggsave("../../../other_repos/GlobalDA_paper/figures/logRchange_breastadenoca.pdf", width = 10, height = 2)

## Now plotting the absolute exposures with these log-ratios

obj_breast
mat_LRchanges = melt(log(normalise_rw(obj_breast$Y[obj_breast$x[,2] == 0,])/normalise_rw(obj_breast$Y[obj_breast$x[,2] == 1,])))
# mat_LRchanges$col = 'NonInf'
# mat_LRchanges$col[is.infinite(mat_LRchanges$value)] = '-Inf'
# mat_LRchanges$value[is.infinite(mat_LRchanges$value)] = min( remove_inf(mat_LRchanges$value))
# mat_LRchanges$col[is.nan(mat_LRchanges$value)] = 'Inf'
# mat_LRchanges$value[is.nan(mat_LRchanges$value)] = max( remove_inf(mat_LRchanges$value))

## two colours: the overall relative abundance in early, and the relative abundance in late
## perhaps sum of relative in early, and sum of relative in late

overall_early = melt(do.call('rbind', replicate(nrow(obj_breast$Y[obj_breast$x[,2] == 0,]), colSums(normalise_rw(obj_breast$Y[obj_breast$x[,2] == 0,])), simplify = F)))
overall_late = melt(do.call('rbind', replicate(nrow(obj_breast$Y[obj_breast$x[,2] == 1,]), colSums(normalise_rw(obj_breast$Y[obj_breast$x[,2] == 1,])), simplify = F)))

mat_LRchanges_all = cbind(overall_early=overall_early,
                          overall_late=overall_late,
      early=melt(obj_breast$Y[obj_breast$x[,2] == 0,]),
      late=melt(obj_breast$Y[obj_breast$x[,2] == 1,]))
mat_LRchanges_all$early.value
mat_LRchanges_all$late.value

library(viridis)
library(viridisLite)
ggplot(mat_LRchanges_all, aes(x=early.value, y=late.value))+
  geom_point(aes(col=log(mat_LRchanges_all$overall_early.value)), size=4)+
  geom_point(aes(col=log(mat_LRchanges_all$overall_late.value)), size=1.5)+
  scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')

## scatterplot in which the x axis is the average exposure of the signature in all samples
##                      the y axis is the fraction of samples with -Inf or Inf in this signature

## fraction of samples with a zero in either group

give_scatterplot_ct = function(ct){
  ct_obj = load_PCAWG(ct, "signatures", simulation=F)
  if(!is.na(ct_obj)){
    fraction_zeros = colSums((ct_obj$Y[ct_obj$x[,2] == 0,] == 0) | (ct_obj$Y[ct_obj$x[,2] == 1,] == 0))/nrow(ct_obj$Y)*2
    
    ## average exposure for nonzero samples
    average_exposure = apply(rbind(normalise_rw(ct_obj$Y[ct_obj$x[,2] == 0,]),
          normalise_rw(ct_obj$Y[ct_obj$x[,2] == 1,])), 2, function(i) mean(i[i != 0]))
    
    return(cbind.data.frame(fraction_zeros, average_exposure, ct, sig=names(fraction_zeros)))
  }
}

all_scatterplot = do.call('rbind', lapply(sort(unique(samples_files$X1)), give_scatterplot_ct))

ggplot(all_scatterplot, aes(x=fraction_zeros, y=average_exposure))+geom_point()
ggsave("../../results/exploratory/scatterplot_sig_exposure_and_zeros.pdf", height = 5, width = 5)

ggplot(all_scatterplot, aes(x=fraction_zeros, y=average_exposure))+geom_point()+facet_wrap(.~sig)
ggsave("../../results/exploratory/scatterplot_sig_exposure_and_zeros_sigfacet.pdf", width = 8, height = 8)

ggplot(all_scatterplot, aes(x=fraction_zeros, y=average_exposure))+geom_point()+facet_wrap(.~ct)
ggsave("../../results/exploratory/scatterplot_sig_exposure_and_zeros_sigCT.pdf", width = 8, height = 8)

pdf("../../results/exploratory/scatterplot_sig_exposure_and_zeros_sigCT_2.pdf")
for(ct_it in unique(all_scatterplot$ct)){
  print(ggplot(all_scatterplot %>% filter(ct == ct_it), aes(x=fraction_zeros, y=average_exposure, label=sig))+
          geom_point()+geom_label()+ggtitle(ct_it))
}
dev.off()

##---------
give_all_but_first_k = function(i,k){i[(k+1):length(i)]}
results_M = (wrapper_run_TMB("Breast-AdenoCA", "signatures", model = "fullRE_M"))
new_obj_subset = give_subset_sigs(readRDS(file = paste0("../../data/roo/", "Breast-AdenoCA", '_', "signatures", "_ROO.RDS" )),
                 sigs_to_remove = give_all_but_first_k(colnames(load_PCAWG("Breast-AdenoCA", "signatures", simulation=F)$Y),4))
new_obj_subset = give_subset_sigs(readRDS(file = paste0("../../data/roo/", "Breast-AdenoCA", '_', "signatures", "_ROO.RDS" )),
                                  sigs_to_remove = c("SBS5", "SBS8", "SBS9", "SBS17a", "SBS17b", "SBS18", "SBS37", "SBS39", "SBS41"))
saveRDS(new_obj_subset,
        file = paste0("../../data/roo/", "Breast-AdenoCA", '-subset_', "signatures", "_ROO.RDS"))
results_DM_subset = (wrapper_run_TMB("Breast-AdenoCA-subset", "signatures", model = "fullRE_DM"))
results_DM_subset
results_DM_subset = (wrapper_run_TMB("Breast-AdenoCA-subset", "signatures", model = "diagRE_DM"))
results_DM_subset

run_wrap = function(){
  results_DM_subset = (wrapper_run_TMB("Breast-AdenoCA-subset", "signatures", model = "fullRE_DM"))
  results_DM_subset
}
replicas = replicate(20, run_wrap(), simplify = F)

give_barplot("Breast-AdenoCA-subset", "signatures", F)

ct = "Lung-SCC"
give_barplot(ct, typedata, F)
results_DM_sparse = (wrapper_run_TMB(ct, typedata, model = "sparseRE_DM")) ## has not converged yet
results_diagRE_M = (wrapper_run_TMB(ct, typedata, model = "diagRE_M")) ## did not converge at first
results_diagRE_DM = (wrapper_run_TMB(ct, typedata, model = "diagRE_DM"))
results_diagRE_DM
results_DM_sparse
results_diagRE_M
length(python_like_select_name(results_DM_sparse$par.fixed, "beta"))/2 ## number of signatures
## crazy number of signatures

pdf("../../results/exploratory/barplot_pcawg.pdf")
for(ct_it in sort(unique(samples_files[,1]))){
  if(!grepl('subset', ct_it)){
    for(type_it in sort(unique(samples_files[,2]))[c(1,3)]){
      give_barplot(ct_it, type_it, F, title=paste0(ct_it, ' ', type_it))
    }
  }
}
dev.off()

#------