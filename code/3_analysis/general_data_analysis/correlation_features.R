
# setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
library(dplyr)
library("pryr")

fles_in = list.files("../data/roo/", full.names=TRUE)
roo_files = sapply(fles_in, readRDS)

fles_in_mod = basename(fles_in)
fles_in_mod = data.frame(do.call('rbind.data.frame', lapply(fles_in_mod, function(i) strsplit(i, '_')[[1]][1:2])))
colnames(fles_in_mod) = c('CT', 'Feature')
feature_ct = c('nucleotidesubstitution1', 'signatures')

unique(fles_in_mod)

names(roo_files)

append = c()
rm_ct = list()

for(ct in fles_in_mod$CT){
  for(feature_ct_it in feature_ct){
    roo = roo_files[[paste0("../data/roo//", ct, '_', feature_ct_it, "_ROO.RDS")]]
    if(!is.na(roo)){
      flename = paste0("../results/reports/figure_roo_correlations/", ct, '_', feature_ct_it, '.png')
      pairs_plt = apply(do.call('rbind', c(slot(roo, "count_matrices_all"))), 2, as.numeric)
      pairs_plt = normalise_rw(pairs_plt)
      png(flename, width = 3, height = 3, units = "in", res = 300)
      if(dim(pairs_plt)[2] < 10){
        pairs(pairs_plt, pch=19, cex=0.4, main=ct)
      }else{
        plot(0, 0, main=ct)
      }
      dev.off()
      if(feature_ct_it == "signatures"){
        flename = paste0("../results/reports/figure_roo_correlations/", ct, '_', feature_ct_it, '_active.png')
        if(length(slot(roo, "count_matrices_active")[[1]]) == 0){
          png(flename, width = 5, height = 5, units = "in", res = 100)
            plot(0, 0, main=ct)
          dev.off()
        }else{
            pairs_plt_active = apply(do.call('rbind', slot(roo, "count_matrices_active")), 2, as.numeric)
            pairs_plt_active = normalise_rw(pairs_plt_active)
            png(flename, width = 3, height = 3, units = "in", res = 100)
            if(dim(pairs_plt)[2] < 10){
              pairs(pairs_plt, pch=19, cex=0.4, main=ct)
            }else{
              plot(0, 0, main=ct)
            }
            dev.off()
        }
      }
      append = c(append, paste0("\\begin{minipage}\\includegraphics[width=3in]{",flename, '}\\end{minipage}' ))
      rm(pairs_plt)
    }else{
      rm_ct = c(rm_ct, ct)
    }
  }
}

names(roo_files)
paste0("../data/roo//", ct, '_', feature_ct_it, "_ROO.RDS")

mem_used()
