##' breast metastases from Caldas (sent by Stephen-John Sammut),
##' paper  The Genomic and Immune Landscapes of Lethal Metastatic Breast Cancer 

setwd("/Users/morril01/Documents/PhD/CDA_in_Cancer/")

data_metastases = read.table("data/brca_metastases/mutations.maf", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

dim(data_metastases)
samples = unique(data_metastases$sample)
patients = unique(sapply(data_metastases$sample, function(i) strsplit(i, '-')[[1]][1]))
data_metastases$patient = sapply(data_metastases$sample, function(i) strsplit(i, '-')[[1]][1])
## There should be ten patients
length(patients)
length(samples)

per_patient = lapply(patients, function(p) subset(data_metastases, data_metastases$patient == p))
names(per_patient) = patients

per_patient$`288`

per_sample = lapply(per_patient, function(p){
  .x = lapply(unique(p$sample), function(s) subset(p, p$sample == s))
  names(.x) = unique(p$sample)
  .x
  })
names(per_sample) = patients
View(per_sample$`288`$`288-005`)
