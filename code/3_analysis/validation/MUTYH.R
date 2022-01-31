read_annotation

source("../../2_inference_TMB/helper_TMB.R")

read_annotation <- read_annotation[sapply(read_annotation, typeof) != "character"]
read_annotation[[1]]

GOI <- lapply(read_annotation, function(i) i[i$gene_symbol %in% c('MUTYH'),])
GOI_muts_NS <- lapply(GOI, function(i) i[i$CONSEQUENCE == "nonsynonymous",])

num_muts_MUTYH <- sapply(GOI_muts_NS, length)
max(num_muts_MUTYH)

GOI_muts_NS[[which.max(num_muts_MUTYH)]]

# get_ct <- 
samples <- gsub("_mutation_ccf.txt", "", basename(names(GOI)))
samples
metadata <- read.table("../../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", head=T)
ct_relevant <- unique(metadata[metadata$samplename %in% samples,]$histology_detailed)


muts <- lapply(ct_relevant, function(ct) load_PCAWG(ct = ct, typedata = "signaturesPCAWG",
                                                    path_to_data = "../../../data/",
                                                    load_all_sigs = T))
mutsMSE <- lapply(ct_relevant, function(ct) try(load_PCAWG(ct = ct, typedata = "signaturesMSE",
                                                    path_to_data = "../../../data/",
                                                    load_all_sigs = T)))
names(muts) <- names(mutsMSE) <- ct_relevant


remove_na <- function(i) i[!is.na(i)]

muts_to_df <- function(muts){
  muts <- lapply(muts, function(i) i['Y'])
  muts <- lapply(muts, function(i) i[[1]])
  mutsrbind <- do.call('rbind', muts)
  mutsrbind_mut <- mutsrbind[match(samples, rownames(mutsrbind)),]
  mutsrbind_nomut <- mutsrbind[-remove_na(match(samples, rownames(mutsrbind))),]
  return(list(mutsrbind_mut=mutsrbind_mut, mutsrbind_nomut=mutsrbind_nomut))
}

res_muts <- muts_to_df(muts)
mutsrbind_mut <- res_muts[['mutsrbind_mut']]
mutsrbind_nomut <- res_muts[['mutsrbind_nomut']]

res_mutsMSE <- muts_to_df(mutsMSE)
mutsrbind_mutMSE <- res_mutsMSE[['mutsrbind_mut']]
mutsrbind_nomutMSE <- res_mutsMSE[['mutsrbind_nomut']]


boxplot(c(mut=list(mutsrbind_mut[,'SBS36']),
          nomut=list(mutsrbind_nomut[,'SBS36'])))
boxplot(c(mut=list(log(mutsrbind_mut[,'SBS36'])),
          nomut=list(log(mutsrbind_nomut[,'SBS36']))))

boxplot(c(mut=list(mutsrbind_mut[,'SBS8']),
          nomut=list(mutsrbind_nomut[,'SBS8'])))
boxplot(c(mut=list(log(mutsrbind_mut[,'SBS8'])),
          nomut=list(log(mutsrbind_nomut[,'SBS8']))))

boxplot(c(mut=list(log( (1+mutsrbind_mut[,'SBS8'])/(1+mutsrbind_mut[,'SBS36']))),
          nomut=list(log((1+mutsrbind_nomut[,'SBS8'])/(1+mutsrbind_nomut[,'SBS36'])))))



boxplot(c(mut=list(as.numeric(mutsrbind_mutMSE[,'SBS36'])),
          nomut=list(as.numeric(mutsrbind_nomutMSE[,'SBS36']))))
boxplot(c(mut=list(log(as.numeric(mutsrbind_mutMSE[,'SBS36']))),
          nomut=list(log(as.numeric(mutsrbind_nomutMSE[,'SBS36'])))))

boxplot(c(mut=list(as.numeric(mutsrbind_mutMSE[,'SBS8'])),
          nomut=list(as.numeric(mutsrbind_nomutMSE[,'SBS8']))))
boxplot(c(mut=list(log(as.numeric(mutsrbind_mutMSE[,'SBS8']))),
          nomut=list(log(as.numeric(mutsrbind_nomutMSE[,'SBS8'])))))

