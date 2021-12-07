flder = "../data/assessing_models_simulation/inference_results/TMB/nlminb/"

a = list.files(flder, full.names = T); a = a[grepl('multiple_', a)]; x = sapply(a, readRDS)


table(unlist(sapply(x, `[`, 'pdHess')))

## delete files in which there was no convergence
#for(i in names(x)[sapply(x, `[`, 'pdHess') != TRUE]){system(paste0('rm ', i, '\n'))}

