## For the samples for which there is a clear clonal-subclonal structure, split them into clonal vs subclonal

debug_bool=FALSE

if(debug_bool){
  rm(list = ls())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  setwd("../")
  Sys.setenv(LANG='en')
  debug_bool=TRUE
}

suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
source("1_create_ROO/helper_1_create_ROO.R")

if(debug_bool){
  opt=list();
  opt$input_files='../data/restricted/pcawg/pcawg_restricted_snv_counts/00b9d0e6-69dc-4345-bffd-ce32880c8eef ../data/restricted/pcawg/pcawg_restricted_snv_counts/02917220-6a7a-46a1-8656-907e96bef88e ../data/restricted/pcawg/pcawg_restricted_snv_counts/03ad38a6-0902-4aaa-84a3-91ea88fa9883 ../data/restricted/pcawg/pcawg_restricted_snv_counts/05616329-e7ba-4efd-87b1-d79cd0f7af3d ../data/restricted/pcawg/pcawg_restricted_snv_counts/068f4f69-d2fe-4f25-912e-ca7d4623efb6 ../data/restricted/pcawg/pcawg_restricted_snv_counts/0e7f46ca-6f5c-4538-b6d6-00af65d57fcf ebe0ed67-2d3f-45cd-8f9b-4912595b16a0 ../data/restricted/pcawg/pcawg_restricted_snv_counts/fa676301-902f-473f-8313-5bff34ae549a ../data/restricted/pcawg/pcawg_restricted_snv_counts/fce8d8c6-f2a0-43a8-9a7a-b9c519a6686c'
  opt$cancer_type=" Lymph-BNHL";
  opt$feature_type="nucleotidesubstitution1"; opt$output="../data/roo/Lymph-BNHL_nucleotidesubstitution1_ROO.RDS"
}else{
  option_list = list(
    make_option(c("--cancer_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--feature_type"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--output"), type="character", default=NA, 
                help="", metavar="character"),
    make_option(c("--input_files"), type="character", default=NA, 
                help="", metavar="character")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
}

type_features_vec = c('nucleotidesubstitution1', 'nucleotidesubstitution3', 'signatures')


if(! (opt$feature_type %in% type_features_vec)){stop('Incorrect value for <feature_type>')}

flder_clonal = "../data/restricted/pcawg/consensus_subclonal_reconstruction_20170325/"
flder_vcf = "../data/restricted/pcawg/pcawg_restricted_snv_counts/"
fles_vcf = list.files(flder_vcf)

fles_cancer_type = paste0(basename(strsplit(opt$input_files, " ")[[1]]), "_cluster_assignments.txt.gz")

if(length(fles_cancer_type) == 0){
  stop()
}

objects_file = mclapply(fles_cancer_type, function(name_file){
  ## first, read in the clonal deconvolution files from PCAWG
  ## then, create the ROO objects with the existing function createRDS_ROOSigs 
  
  gunzip(paste0(flder_clonal,name_file ))                              ## unzip
  filename_unzipped = paste0(flder_clonal, gsub(".gz", "", name_file)) 
  clonal_file = read.table(filename_unzipped, header = TRUE)
  gzip(filename_unzipped)                                              ## zip back again
  
  raw_name = gsub('_cluster_assignments.txt', '', basename(filename_unzipped))
  VAF_CCF_filename = fles_vcf[grep(raw_name, fles_vcf)]
  if(length(VAF_CCF_filename) == 0){
    ## no VAF file found.
    cat('File ', raw_name, '(cancer type=',opt$cancer_type, ') does not have a VCF file')
    return(list(NA,NA, "No VAF file"))
  }
  name_file_vcf = paste0(flder_vcf, VAF_CCF_filename)
  VAF_CCF_file = read.table(name_file_vcf, sep = "\t", header = TRUE)
  
  ## merge
  merged = merge(clonal_file, VAF_CCF_file, by.x=c("chromosome", "position"), by.y=c("chromosome", "position"))
  merged[,'bool_group_1'] = merged$cluster_1 > 0.5 ## boolean for the two groups
  merged = merged[!(is.na(merged[,'bool_group_1'])),]
  ## now split with the posterior of cluster_1 vs cluster_2, and get in the same format has before
  
  if(length(table(merged$bool_group_1)) == 1){
    # if(length(fles_cancer_type) == 1){
    #   saveRDS(NA, file = opt$output) ## there is only one sample
    #   quit()
    # }else{
      return(list(NA,NA, "Empty group #1"))
    # }
  }
  
  ## Read in the cancer types
  # opt$cancer_type = pcawg_data[grep(gsub('.cluster_assignments.txt.gz', '', name_file), pcawg_data$File.Name),'Project']
  
  ###################################################################################################
  #### Split mutations into groups (e.g. clonal and subclonal), and save as RDS (or return file) ####
  ###################################################################################################
  # cat('Getting RDS object for file ', name_file)
  .tmp_roo_obj = createRDS_ROOSigs_object(pre_path="1_create_ROO/",
                    vcf_path="../data/restricted/pcawg/pcawg_restricted_snv/",
                    ccf_threshold=NA,
                    type_features=opt$feature_type,
                    outfolder=basename(opt$output),
                    active_sigs_version="active_signatures_transposed_clonesig2",
                    cancer_type=opt$cancer_type,
                    in_dataframe=merged)
  # if(is.null(.tmp_roo_obj)){
  #   ## There are no signatures, because of a problem
  #   .tmp_roo_obj = list(list(list(), list()), list(list(), list()))
  # }
  
  if(opt$feature_type == "signatures"){
  #   if(any(sapply(.tmp_roo_obj[[1]], length) == 0)){
  #     return(list(NA,NA, "Empty group #2"))
  #   }
    rds_object = .tmp_roo_obj[[1]] ## this is a list of 2 of 2 lists, subsetted
    rds_object_active_sigs = .tmp_roo_obj[[2]] ## this is a list of 2 of 2 lists, subsetted
  }else{
  #   ## if one group only has zeros, don't return
  #   if(any(sapply(.tmp_roo_obj, sum) == 0)){
  #     return(list(NA,NA, "Empty group #2"))
  #   }
    rds_object = .tmp_roo_obj ## this is a list of 1
    rds_object_active_sigs = list(list(), list())
  }
  
  return(list(rds_object, rds_object_active_sigs, "successful"))
  # stepi = stepi + 1
})

successful_bool = ( sapply(objects_file, function(i) i[[3]]) == "successful")
rds_object = lapply(objects_file, function(i) i[[1]])
rds_object_active_sigs = lapply(objects_file, function(i) i[[2]])

names(rds_object) = names(rds_object_active_sigs) = gsub("_cluster_assignments.txt.gz", "", fles_cancer_type);
rds_object = rds_object[successful_bool]
rds_object_active_sigs = rds_object_active_sigs[successful_bool]


#' in the case that we don't have active signatures - either because we are looking at mutation categories as
#' opposed to signatures, or because active signatures are not defined, rds_object_active_sigs should be of
#' the form list(list(), list())

## for all the files in a cancer type
###########################################
#### Creating the ROO objects from RDS ####
###########################################

## Now we want to put together all the objects of the relevant cancer type
if(!(length(rds_object) == 0)){
  all_objs_activesigs = rds_object_active_sigs
  if(opt$feature_type=="signatures"){
    all_objs_activesigs = rds_object_active_sigs

    ## if the cancer type does not has active signatures, change the name of the files
    if(all(sapply( all_objs_activesigs, function(i) (length(i[[1]]) == 0) & (length(i[[2]]) == 0) ))){
      names(all_objs_activesigs) = NULL
    }
  }
  robject = createROO_ROOSigs_object(type_features=opt$feature_type,
                                 all_objs=rds_object,
                                 all_objs_activesigs=all_objs_activesigs,
                                 save_bool=FALSE,
                                 outfiles=names(rds_object),
                                 pre_path="1_create_ROO/",
                                 cancer_type_given = rep(opt$cancer_type, length(rds_object)))
  saveRDS(robject, file = opt$output)
}else{
  saveRDS(NA, file = opt$output)
}

