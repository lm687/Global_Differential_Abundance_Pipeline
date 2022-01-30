##' read ccf files and see if mutations of certain pathways are mutated predominantly in clonal
##' or subclonal
##' 

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

flder <- "../../../data/restricted/pcawg/consensus_subclonal_reconstruction_mutccf_20170325/"

library(VariantAnnotation)
library(annotate)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

source("helper_load_vcf_and_clonal.R")

change_column_names <- function(i){
  colnames(i)[colnames(i) == "startpos"] <- "start"
  colnames(i)[colnames(i) == "endpos"] <- "end"
  i
}


add_start_end <- function(i){
  colnames(i)[colnames(i) == "position"] <- "start"
  colnames(i)[colnames(i) == "position.1"] <- "end"
  i
}

input_files <- list.files(flder, full.names = T)

already_computed <- gsub('.RDS', '', list.files("../../../data/restricted/pcawg/var_annotation_with_ccf"))
bool_already_computed <- sapply(input_files, function(i) any(sapply(already_computed, function(j) grepl(j, i))))
table(bool_already_computed)

for(file_name in input_files[!bool_already_computed]){
  basename(file_name)
  # 0009b464-b376-4fbc-8a56-da538269a02f
  # 0009b464-b376-4fbc-8a56-da538269a02f
  fle_vcf_path <- paste0("../../../data/restricted/pcawg/pcawg_restricted_snv/", gsub("_mutation_ccf.txt", "", basename(file_name)), ".consensus.20160830.somatic.snv_mnv.vcf.gz")
  if(file.exists(fle_vcf_path)){
    vcf_file <- fle_vcf_path
  }else{
    vcf_file <- gsub(".gz", "", fle_vcf_path) ## try to read if it's unzipped
  }
  if(!(file.exists(fle_vcf_path))){
    next
  }
  # file2 <- read.table( paste0("../../../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/", gsub("_mutation_ccf.txt", "", basename(file_name)), "_mutation_timing.txt.gz"), se='\t', header=T)

  .matched_files <- give_unified_vcf_and_ccf(vcf_file, file_name)
  subclonality_df <- .matched_files$subclonality_df
  filemut <- .matched_files$vcf
  
  # loc_all <- locateVariants(.query, txdb, AllVariants())

  .query <- as(data.frame(add_start_end(cbind.data.frame(subclonality_df[,c('chromosome', 'position', 'position')], segVal=NA))), "GRanges")
  coding_predict <- predictCoding(query = .query, subject = txdb, seqSource = Hsapiens,
                                  varAllele = DNAStringSet(vcfR::getALT(filemut)))

  if(length(coding_predict) > 0){
    .symbols <- getSYMBOL(c(coding_predict$GENEID[!is.na(coding_predict$GENEID)]), data='org.Hs.eg')
    coding_predict$gene_symbol[!is.na(coding_predict$GENEID)] = .symbols
    coding_predict$ccf <- subclonality_df$ccf[coding_predict$QUERYID]
  }
  # View(as(coding_predict, 'data.frame'))
  saveRDS(coding_predict, file = gsub(".txt", ".RDS", gsub(".consensus_subclonal_reconstruction_mutccf_20170325/", "/var_annotation_with_ccf", file_name)))
  # unique(coding_predict$gene_symbol)
  
}

# 
# 
# library(GenomicRanges)
# library(ensembldb)
# library(parallel)
# 
# gtf.file <- file.path("../../../../other_repos/Vias_Brenton/RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/Homo_sapiens.GRCh37.87.gtf.gz")
# sqlite_file <- 'Homo_sapiens.GRCh37.87.sqlite'
# path_sq <- "../../../../other_repos/Vias_Brenton/RNASeq_and_CN/20191218_ViasM_BJ_orgaBrs/Data/"
# sqlite_path <- file.path(path_sq, sqlite_file)
# 
# if(!file.exists(sqlite_path)) {
#   ## generate the SQLite database file
#   ensembldb::ensDbFromGtf(gtf=gtf.file, path = path_sq, outfile=sqlite_file)
# }
# EnsDb.Hsapiens.v87 <- ensembldb::EnsDb(sqlite_path)
# 
# # Genes, used to annotated the TPM matrix to send to Maria
# ag <- ensembldb::genes(EnsDb.Hsapiens.v87, filter=list(AnnotationFilter::GeneBiotypeFilter('protein_coding')), return.type="DataFrame") 
# ag
# 
# chromosomes_metadata <- ensembldb::seqinfo(EnsDb.Hsapiens.v87)
# chromlens <- (cbind.data.frame(Chrom=seqnames(chromosomes_metadata),
#                                Length=seqlengths(chromosomes_metadata)))
# head(chromlens)
# 
# ag_subsetchrom <- ag[!(ag$seq_name %in% c("MT", "X", "Y")) & !grepl("GL",ag$seq_name),]
# unique(ag_subsetchrom$seq_name)
# gr_genes = as(paste0(ag_subsetchrom$seq_name, ':', ag_subsetchrom$gene_seq_start, '-', ag_subsetchrom$gene_seq_end), "GRanges")
# gr_genes
# 
# 
# levels(seqnames(gr_genes))
# seqnames(gr_genes) <- droplevels(seqnames(gr_genes))
# levels(seqnames(gr_genes))
# unique(seqnames(gr_genes))
# 
# give_CN_per_gene <- function(segment_arg){
#   ## in this case GR_bins is gr_genes
#   
#   GR_bins <- gr_genes
#   gr_CN = as(data.frame(segment_arg), "GRanges")
#   GR_bins <- trim(GR_bins)
#   
#   gr_CN <- trim(gr_CN)
#   
#   values(gr_CN) = segment_arg[,'segVal']
#   
#   ## do per sample
#   gr_CN_org = gr_CN
#   seqlevels(gr_CN) ## necessary for mcolAsRleList
#   seqlengths(gr_CN) = as.numeric(as.vector(chromlens$Length[match(seqlevels(gr_CN), gsub("chr", "", chromlens$Chrom))])) ## necessary for mcolAsRleList
#   ag_subsetchrom_subset <- ag_subsetchrom[which(seqnames(GR_bins) %in% seqlevels(gr_CN)),]
#   GR_bins <- GR_bins[seqnames(GR_bins) %in% seqlevels(gr_CN)]
#   seqlevels(GR_bins) = seqlevels(gr_CN) ## necessary for mcolAsRleList
#   seqlengths(GR_bins) = seqlengths(gr_CN) ## necessary for mcolAsRleList
#   
#   
#   full_GR <- c(GR_bins, gr_CN_org)
#   disjoint_gr <- GenomicRanges::disjoin(full_GR, with.revmap=TRUE, ignore.strand=TRUE)
#   disjoint_gr <- trim(disjoint_gr)
#   
#   disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
#   
#   ## keep only those in which there is an overlap with a gene, i.e. revmap has to include the idxs 1:length(gr_genes)
#   disjoint_gr = disjoint_gr[disjoint_gr_revmap_first %in% 1:length(GR_bins),]
#   disjoint_gr_revmap_first = sapply(disjoint_gr$revmap, function(i) i[1])
#   disjoint_gr_revmap_second = sapply(disjoint_gr$revmap, function(i) i[2])
#   table(is.na(disjoint_gr_revmap_second))
#   
#   ## remove sections which only contain non-genes
#   disjoint_gr[is.na(disjoint_gr_revmap_second),]
#   
#   ## get the idx of the copy number segments
#   idx_gr = disjoint_gr_revmap_second-length(GR_bins)
#   is.na(idx_gr)
#   
#   idx_gr[idx_gr <= 0 ] = NA
#   
#   disjoint_gr$CN_val = gr_CN_org$X[idx_gr]
#   disjoint_gr$CN_val[is.na(disjoint_gr$CN_val)] = 2
#   
#   disjoint_gr$revmap = disjoint_gr_revmap_first
#   seqlevels(disjoint_gr) = as.character(unique(seqnames(disjoint_gr)))
#   seqlengths(disjoint_gr) = as.numeric(as.vector(chromlens$Length[match(seqlevels(disjoint_gr), gsub("chr", "", chromlens$Chrom))]))
#   seqlevels(GR_bins) = seqlevels(disjoint_gr)
#   seqlengths(GR_bins) = seqlengths(disjoint_gr)
#   
#   CN_RleList = mcolAsRleList(disjoint_gr, "CN_val") ## that takes forever even when length(disjoint_gr) is a VERY low value
#   averageCN = binnedAverage(bins = GR_bins, numvar = CN_RleList, varname = "CN_bin_averaged")
#   
#   return(list(averageCN=averageCN, GR_bins=ag_subsetchrom_subset))
# }
# 
# ## get segments from the unified version created in copy_number_analysis_organoids.Rmd
# segments_Britroc <- readRDS("../../../../other_repos/Vias_Brenton/copy_number_analysis_organoids/data/clean_segtables/segtables_BriTROC_absolute_copynumber.RDS")
# segments_Britroc <- segments_Britroc[1:2]
# 
# 
# names_not_run
# 
# head(segments_Britroc[[1]])
# 
# segments_Britroc <- list(example=add_start_end(cbind.data.frame(head(a[,c('chromosome', 'position', 'position')]), segVal=NA)))
# sample_name_it <- names(segments_Britroc)
# 
# CN_averages = lapply(names(segments_Britroc), function(sample_name_it){
#   CN_averages <- give_CN_per_gene(segment_arg = change_column_names(segments_Britroc[[sample_name_it]]))
#   averaged_CN_df = cbind.data.frame(gene_name=CN_averages$GR_bins$symbol,
#                                     CN=CN_averages$averageCN$CN_bin_averaged)
#   saveRDS(averaged_CN_df, paste0("../output/output_GRCh37/all_CN_states_per_gene/", sample_name_it, ".RDS"))
# })

