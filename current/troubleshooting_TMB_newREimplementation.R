### Troubleshooting inference with TMB

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

ct = "Kidney-RCC.clearcell"
ct = "Thy-AdenoCA"
ct = "SoftTissue-Liposarc"
ct = "Uterus-AdenoCA"
typedata =  "signatures" #samples_files[1,2]
samples_files = data.frame(do.call('rbind', sapply(gsub("_ROO.RDS", "", list.files("../../data/roo/")),
                                                   strsplit, split = "_")))

TMB::compile("mm_multinomial/fullRE_dirichletmultinomial_single_lambda_REv2.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_dirichletmultinomial_single_lambda_REv2"))
TMB::compile("mm_multinomial/fullRE_ME_multinomial_REv2.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial_REv2"))
TMB::compile("mm_multinomial/fullRE_ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_multinomial"))
TMB::compile("mm_multinomial/fullRE_ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/fullRE_ME_dirichletmultinomial"))

## it's 5-dimensional
fivedimobject = load_PCAWG(ct = ct, typedata = typedata)
# fivedimobject = give_subset_sigs_TMBobj(fivedimobject, sigs_to_remove = c('SBS3', 'SBS24', 'SBS26', 'SBS22', 'SBS29', 'SBS40', 'SBS41', 'SBS13'))
fivedimobject = give_subset_sigs_TMBobj(fivedimobject,
  sigs_to_remove = names(sort(colSums(fivedimobject$Y), d=T)[-(1:5)]))

stopifnot(ncol(fivedimobject$Y) == 5)

results_DMsl = (wrapper_run_TMB(object = fivedimobject, model = 'fullRE_dirichletmultinomial_singlelambda_REv2'))
results_DMsl
results_DMoldRE = (wrapper_run_TMB(object = fivedimobject, model = 'fullRE_DM'))
results_DMoldRE
results_M = (wrapper_run_TMB(object = fivedimobject, model = 'fullRE_multinomial_REv2'))
results_M
results_MoldRE = (wrapper_run_TMB(object = fivedimobject, model = 'fullRE_M'))
results_MoldRE

cbind_est <- cbind.data.frame(results_DMsl=results_DMsl$par.fixed[1:18],
                              results_DMoldRE=results_DMoldRE$par.fixed[1:18],
                              results_M=results_M$par.fixed,
                              results_MoldRE=results_MoldRE$par.fixed)
# pairs(cbind_est)

plts= replicate(n = ncol(cbind_est)^2, list())
for(i in 1:ncol(cbind_est)){
  for(j in 1:ncol(cbind_est)){
    if(i != j){
      plt_it <- ggplot(cbind.data.frame(x=cbind_est[,i], y=cbind_est[,j], col=names(results_DMsl$par.fixed[1:18])),
             aes(x=x,y=y, col=col))+geom_point()+labs(x=colnames(cbind_est)[i], y=colnames(cbind_est)[j])+
        geom_abline(slope = 1, intercept = 0, lty='dashed')
    }else{
      plt_it <- ggplot()
    }
    plts[[(i-1)*4+j]] = plt_it
  }
}

do.call("grid.arrange", list(grobs=plts, ncol=4))

