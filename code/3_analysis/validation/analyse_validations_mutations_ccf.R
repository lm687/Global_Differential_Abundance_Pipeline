setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(reshape2)
library(ggplot2)
library(dplyr)

flder <- "../../../data/restricted/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/"

metadata <- read.table("../../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", head=T)
enough_samples = read.table("../../../data/restricted/pcawg/CT_sufficient_samples.txt", comment.char='#')[,1]

input_files <- list.files(flder, full.names = T)

file_name <- input_files
a <- sapply(file_name, function(i) gsub(".txt", ".RDS", gsub(".consensus_subclonal_reconstruction_mutccf_20170325/", "/var_annotation_with_ccf", i)))

read_annotation <- lapply(a, function(j) try(readRDS(j)))
read_annotation <- read_annotation[sapply(read_annotation, typeof) != "character"]

# kegg_pathway_name <- 'KEGG HRD'
# kegg_pathway_name <- 'KEGG MMR'
# kegg_pathway_name <- 'KEGG BER'
# kegg_pathway_name <- 'KEGG NER'
# kegg_pathway_name <- 'KEGG NHEJ'
kegg_pathway_name <- 'TP53'
single_gene <- F
if(kegg_pathway_name == 'KEGG HRD'){
  kegg <- read.table("../../../data/other/KEGG/KEGG_HOMOLOGOUS_RECOMBINATION.txt", sep = "\t")
}else if(kegg_pathway_name == 'KEGG MMR'){
  kegg <- read.table("../../../data/other/KEGG/KEGG_MISMATCH_REPAIR.txt", sep = "\t")
}else if(kegg_pathway_name == 'KEGG BER'){
  kegg <- read.table("../../../data/other/KEGG/KEGG_BASE_EXCISION_REPAIR.txt", sep = "\t")
}else if(kegg_pathway_name == 'KEGG NER'){
  kegg <- read.table("../../../data/other/KEGG/KEGG_NUCLEOTIDE_EXCISION_REPAIR.txt", sep = "\t")
}else if(kegg_pathway_name == 'KEGG NHEJ'){
  kegg <- read.table("../../../data/other/KEGG/KEGG_NON_HOMOLOGOUS_END_JOINING.txt", sep = "\t")
}else if(kegg_pathway_name == 'TP53'){
  kegg <- 'TP53'
  single_gene <- T
}else{
  
}

if(single_gene){
  kegg <- 'CCNE1'
  kegg_pathway_name <- kegg
  keggname_out <- gsub(" ", "", kegg_pathway_name)
}

if(!single_gene) kegg <- kegg$V1[-c(1:2)] ## first two lines are annotation
GOI <- lapply(read_annotation, function(i) i[i$gene_symbol %in% kegg,])

GOI_len <- sapply(GOI, length)
table(GOI_len)

GOI_muts <- GOI[GOI_len> 0]

GOI_muts_NS <- sapply(GOI_muts, function(i) i[i$CONSEQUENCE == "nonsynonymous",])
GOI_muts_NS_ccf <- lapply(GOI_muts_NS, function(i) data.frame(i)[,c('ccf', 'gene_symbol')])
length(GOI_muts_NS_ccf)
hist(unlist(GOI_muts_NS_ccf))

GOI_muts_NS_ccf_melt <- cbind.data.frame(sample=rep(basename(names(GOI_muts_NS_ccf)), sapply(GOI_muts_NS_ccf, nrow)),
                                          ccf=do.call('rbind', GOI_muts_NS_ccf))
rownames(GOI_muts_NS_ccf_melt) <- NULL
GOI_muts_NS_ccf_melt$ct <- metadata$histology_detailed[match(gsub("_mutation_ccf.txt", "", GOI_muts_NS_ccf_melt$sample),
                                                              metadata$samplename)]
head(GOI_muts_NS_ccf_melt)

metadata

GOI_muts_NS_ccf_melt_dplyr <- GOI_muts_NS_ccf_melt %>% group_by(ct) %>% summarise(median(ccf.ccf))

GOI_muts_NS_ccf_melt <- GOI_muts_NS_ccf_melt[GOI_muts_NS_ccf_melt$ct %in% enough_samples,]

unique(GOI_muts_NS_ccf_melt$ct) %in% GOI_muts_NS_ccf_melt_dplyr$ct
GOI_muts_NS_ccf_melt <- GOI_muts_NS_ccf_melt[!(is.na(GOI_muts_NS_ccf_melt$ct)),]
ggplot(GOI_muts_NS_ccf_melt, aes(x=factor(ct, levels=GOI_muts_NS_ccf_melt_dplyr$ct[order(GOI_muts_NS_ccf_melt_dplyr$`median(ccf.ccf)`)]),
                                  y=ccf.ccf,col=ccf.gene_symbol, group=ct))+geom_boxplot()+geom_jitter(alpha=0.2)+theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='', y= paste0('CCF of non-synonymous ', kegg_pathway_name, ' mutation'), col='Gene')
ggsave(paste0("../../../results/validation/NS_", keggname_out, "_ccf.pdf"), height = 4.5)

GOI_muts_NS_ccf_melt$ccf_lim <- sapply(GOI_muts_NS_ccf_melt$ccf.ccf, function(i) min(i, 1))
ggplot(GOI_muts_NS_ccf_melt, aes(x=factor(ct, levels=GOI_muts_NS_ccf_melt_dplyr$ct[order(GOI_muts_NS_ccf_melt_dplyr$`median(ccf.ccf)`)]),
                                  y=ccf_lim, col=ccf.gene_symbol, group=ct))+geom_boxplot()+geom_jitter(alpha=0.2)+theme_bw()+
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust=1))+
  labs(x='', y= paste0('CCF of non-synonymous ', kegg_pathway_name, ' mutation'), col='Gene')
ggsave(paste0("../../../results/validation/NS_",keggname_out, "_ccf_lim1.pdf"), height = 4.5)

GOI_muts_NS_ccf_melt[is.na(GOI_muts_NS_ccf_melt$ct),]

sort(unique(GOI_muts_NS_ccf_melt_dplyr$ct))
sort(unique(metadata$histology_detailed))

some_gbm <- read_annotation[paste0('../../../data/restricted/pcawg/consensus_subclonal_reconstruction_mutccf_20170325//',
       metadata$samplename[metadata$histology_detailed == 'CNS-GBM'][1:2],
       '_mutation_ccf.txt')]

View(data.frame(some_gbm[[1]][order(some_gbm[[1]]$ccf),]))

'b2ec0fd0-fbcf-4abc-ad80-4ae444e30b55' %in% metadata$samplename
'2df02f2b-9f1c-4249-b3b4-b03079cd97d9' %in% metadata$samplename
