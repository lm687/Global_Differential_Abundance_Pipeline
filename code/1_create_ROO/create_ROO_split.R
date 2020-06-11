## For the samples for which there is a clear clonal-subclonal structure, split them into clonal vs subclonal

debug_bool=TRUE

if(debug_bool){
  rm(list = ls())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  Sys.setenv(LANG='en')
  debug_bool=TRUE
}

suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
source("helper_1_create_ROO.R")

if(debug_bool){

  opt=list(); opt$pcawg_data="../../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt";
  opt$cancer_type="Skin-Melanoma.mucosal";
  opt$feature_type="nucleotidesubstitution3"; opt$output="../../data/roo/Skin-Melanoma.mucosal_nucleotidesubstitution3_ROO.RDS"
}else{
  option_list = list(
    make_option(c("--cancer_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--feature_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--output"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--pcawg_data"), type="character", default=NA, 
                help="", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

type_features_vec = c('nucleotidesubstitution1', 'nucleotidesubstitution3', 'signatures')


if(! (opt$feature_type %in% type_features_vec)){stop('Incorrect value for <feature_type>')}

flder_clonal = "../../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/"
flder_vcf = "../../data/restricted/pcawg/pcawg_restricted_snv_counts/"
fles_vcf = list.files(flder_vcf)

## Run this for each cancer type separately
pcawg_data <- read.table(opt$pcawg_data,
                         stringsAsFactors = FALSE, sep = '\t', header = T)

pcawg_data <- pcawg_data[!duplicated(pcawg_data$samplename),]
all_files = (list.files("~/Documents/PhD/CDA_in_Cancer/data/pcawg/consensus_subclonal_reconstruction_20170325/"))
all_files = subset(all_files, grepl('cluster_assignments', all_files))
pcawg_data_ibj2 = as.vector(sapply(pcawg_data$samplename, function(i) paste0(strsplit(i, '[.]')[[1]][1], collapse="-")))
mrdge = match(gsub("_cluster_assignments.txt.gz", "", all_files), pcawg_data_ibj2)
pcawg_data = pcawg_data[!is.na(mrdge),]
pcawg_data <- pcawg_data[order(pcawg_data$histology_detailed),]

pcawg_data_subset = subset(pcawg_data, pcawg_data$histology_detailed == opt$cancer_type)
pcawg_data_ibj2 = as.vector(sapply(pcawg_data_subset$samplename, function(i) paste0(strsplit(i, '[.]')[[1]][1], collapse="-")))
fles_cancer_type = all_files[match( pcawg_data_ibj2, gsub("_cluster_assignments.txt.gz", "", all_files))]
fles_cancer_type = fles_cancer_type[!is.na(fles_cancer_type)]

rds_object = list()
active_sigs = list()
failed_name_files= c()

if(length(fles_cancer_type) == 0){
  stop()
}


for(name_file in fles_cancer_type){
  ## first, read in the clonal deconvolution files from PCAWG
  ## then, create the ROO objects with the existing function createRDS_ROOSigs 
  
  gunzip(paste0(flder_clonal,name_file ))                              ## unzip
  filename_unzipped = paste0(flder_clonal, gsub(".gz", "", name_file)) 
  clonal_file = read.table(filename_unzipped, header = TRUE)
  gzip(filename_unzipped)                                              ## zip again
  
  raw_name = gsub('_cluster_assignments.txt', '', basename(filename_unzipped))
  VAF_CCF_filename = fles_vcf[grep(raw_name, fles_vcf)]
  if(length(VAF_CCF_filename) == 0){
    ## no VAF file found.
    cat('File ', raw_name, '(cancer type=',opt$cancer_type, ') does not have a VCF file')
    failed_name_files = c(failed_name_files, name_file)
    next
  }
  name_file_vcf = paste0(flder_vcf, VAF_CCF_filename)
  
  VAF_CCF_file = read.table(name_file_vcf, sep = "\t", header = TRUE)
  
  ## merge
  merged = merge(clonal_file, VAF_CCF_file, by.x=c("chromosome", "position"), by.y=c("chromosome", "position"))
  merged[,'bool_group_1'] = merged$cluster_1 > 0.5
  merged = merged[!(is.na(merged[,'bool_group_1'])),]
  ## now split with the posterior of cluster_1 vs cluster_2, and get in the same format has before
  
  if(length(table(merged$bool_group_1)) == 1){
    next ## there is only one group
  }
  
  in_dataframe=merged
  
  ## Read in the cancer types
  # opt$cancer_type = pcawg_data[grep(gsub('.cluster_assignments.txt.gz', '', name_file), pcawg_data$File.Name),'Project']
  
  ###################################################################################################
  #### Split mutations into groups (e.g. clonal and subclonal), and save as RDS (or return file) ####
  ###################################################################################################
  cat('Getting RDS object for file ', name_file)
  .tmp_roo_obj = createRDS_ROOSigs_object(pre_path="",
                    vcf_path="../../data/restricted/pcawg/pcawg_restricted_snv/",
                    ccf_threshold=NA,
                    type_features=opt$feature_type,
                    outfolder=basename(opt$output),
                    active_sigs_version="active_signatures_transposed_clonesig2",
                    cancer_type=opt$cancer_type,
                    in_dataframe=merged)
  if(is.null(.tmp_roo_obj)){
    ## There are no signatures, because of a problem
    .tmp_roo_obj = list(list(list(), list()), list(list(), list()))
  }
  

  if(opt$feature_type == "signatures"){
    if(any(sapply(.tmp_roo_obj[[1]], length) == 0)){
      failed_name_files = c(failed_name_files, name_file)
      next
    }
    rds_object[[name_file]] = .tmp_roo_obj[[1]] ## this is a list of 2 of 2 lists, subsetted
    active_sigs[[name_file]] = .tmp_roo_obj[[2]] ## this is a list of 2 of 2 lists, subsetted
  }else{
    ## if one group only has zeros, don't return
    if(any(sapply(.tmp_roo_obj, sum) == 0)){
      failed_name_files = c(failed_name_files, name_file)
      next
    }
    rds_object[[name_file]] = .tmp_roo_obj[[1]] ## this is a list of 2
  }
  
}

## remove the failed ones??
for(type_features_vec_it in type_features_vec){
  names(rds_object) = gsub("_cluster_assignments.txt.gz", "", names(rds_object))
}

if(opt$feature_type == "signatures"){
  names(active_sigs) = paste0(gsub("_cluster_assignments.txt.gz", "", names(active_sigs)), '_active')
}

## for all the files in a cancer type
###########################################
#### Creating the ROO objects from RDS ####
###########################################

## Now we want to put together all the objects of the relevant cancer type
if(!(length(rds_object) == 0)){
  if(opt$feature_type=="signatures"){
    all_objs_activesigs = active_sigs
    
    ## if the cancer type does not has active signatures, change the name of the files
    if(all(sapply( all_objs_activesigs, function(i) (length(i[[1]]) == 0) & (length(i[[2]]) == 0) ))){
      names(all_objs_activesigs) = NULL
    }
  }else{
    ## Empty list if looking at features
    all_objs_activesigs = list(list(), list())
  }
  robject = createROO_ROOSigs_object(type_features=opt$feature_type,
                                 all_objs=rds_object,
                                 all_objs_activesigs=all_objs_activesigs,
                                 save_bool=FALSE,
                                 outfiles=names(rds_object),
                                 pre_path="",
                                 cancer_type_given = rep(opt$cancer_type, length(rds_object)))
  saveRDS(robject, file = opt$output)
}else{
  saveRDS(NA, file = opt$output)
}



