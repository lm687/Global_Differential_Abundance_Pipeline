flder = "../data/assessing_models_simulation/inference_results/TMB/nlminb/"

a = list.files(flder, full.names = T); a = a[grepl('multiple_GenerationMGnorm_80_180_100_6_0_diagREDM_betaintercept3_betaslope3_sdRE1_dataset', a)]; x = sapply(a, readRDS)


table(unlist(sapply(x, `[`, 'pdHess')))

## delete files in which there was no convergence
#for(i in names(x)[sapply(x, `[`, 'pdHess') != TRUE]){system(paste0('rm ', i, '\n'))}


