rm(list = ls())
setwd("/Users/morril01/Documents/PhD/GlobalDA/code/")
require(reshape2)
require(ggplot2)
require(ggrepel)
folder_objects="../data/roo/"

#----------------------------------------------------------------------------------------#
## Read all the ROO files which contain exposures in two groups
fles_in = list.files(folder_objects, full.names=TRUE)
roo_files = sapply(fles_in, readRDS)
#----------------------------------------------------------------------------------------#

roo_files = roo_files[grep('signatures', names(roo_files))]
names(roo_files) = gsub("_signatures_ROO.RDS", "", basename(names(roo_files)))

get_signatures = function(roo_obj){
  colnames(attr(roo_obj,"count_matrices_active")[[1]])
}

get_signatures(roo_files$`Breast-DCIS`)
sapply(roo_files, get_signatures)


