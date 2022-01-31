give_unified_vcf_and_ccf <- function(vcf_file_path, subclonality_df_path, type_clonal_file='ccf'){
  filemut <- vcfR::read.vcfR(vcf_file_path)
  # bgzip(file = vcf_file, dest = gsub("/pcawg_restricted_snv/", "/pcawg_restricted_snv_bgzip/", vcf_file))
  # filemut <- VariantAnnotation::readVcf(file = "../../../data/restricted/pcawg/pcawg_restricted_snv/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv.vcf.gz", ) ## problem with bgzip
  subclonality_df <- read.table(subclonality_df_path, header = T)
  if(type_clonal_file == 'ccf'){
    subclonality_df <- subclonality_df[subclonality_df$type == 'SNV',] ## no indels
  }else if(type_clonal_file == 'cluster'){
    subclonality_df <- subclonality_df[subclonality_df$mut_type == 'SNV',] ## no indels
  }
  vcf_file_path
  
  .must1 <- paste0(gsub("chr", "", subclonality_df$chromosome), '-', subclonality_df$position)
  .must2 <- paste0( vcfR::getCHROM(filemut), '-',  vcfR::getPOS(filemut))
  table(.must1 %in% .must2)
  table(.must2 %in% .must1)
  .must1[!(.must1 %in% .must2)] ## in clonal file but not in vcf. Vcf is consensus
  .must2[!(.must2 %in% .must1)] ## in vcf but not in clonal file
  
  subclonality_df <- subclonality_df[(.must1 %in% .must2),]
  filemut <- filemut[(.must2 %in% .must1),]
  dim(filemut)
  dim(subclonality_df)
  
  .must1 <- paste0(gsub("chr", "", subclonality_df$chromosome), '-', subclonality_df$position)
  .must2 <- paste0( vcfR::getCHROM(filemut), '-',  vcfR::getPOS(filemut))
  
  ##' in some summarised subclonal reconstruction files only one of the multiple ALT allelles,
  ##' in the same position, is found, whereas in the VCF they all appear
  cat('Mutation without VAF found in vcf but not in annotated clonal file\n')
  subclonality_df[subclonality_df$position == ('50992099'),]
  filemut[vcfR::getPOS(filemut) == '50992099',]@fix
  # file2[file2$position == 50992099,]
  
  sort(table(.must1), decreasing = T)[1:3]
  sort(table(.must2), decreasing = T)[1:3]
  
  .must2[c(which(duplicated(.must2)), which(duplicated(.must2, fromLast=T)))]
  filemut[vcfR::getPOS(filemut) == '50992099',]@fix
  ## select the VCF item that has a VAF value
  idx_missing_from_must1 <- unique(c( which(duplicated(.must2, fromLast=T)), which(duplicated(.must2))))
  
  if(length(idx_missing_from_must1)>0)  info_same_position <- vcfR::getINFO(filemut)[idx_missing_from_must1]
  if(length(idx_missing_from_must1)>0)  filemut <- filemut[-idx_missing_from_must1[!sapply(info_same_position, function(i) grepl('VAF=', i))],]
  
  dim(filemut)
  dim(subclonality_df)
  
  subclonality_df$chromosome <- paste0('chr', subclonality_df$chromosome)
  
  return(list(subclonality_df=subclonality_df, vcf=filemut))
  
}
