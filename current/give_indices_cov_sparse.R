## for each cancer type under consideration, write


sigs_to_keep_cov <- read.table("../../current/subset_sigs_sparse_cov.txt", sep = "\t")
sigs_to_keep_cov <- sigs_to_keep_cov[!(sigs_to_keep_cov$V3 == "?"),]

idx_keep = sapply(sigs_to_keep_cov$V2, function(ct){
  obj_it = load_PCAWG(ct, typedata = "signatures")
  sigs_to_keep_cov_it <- strsplit(sigs_to_keep_cov[sigs_to_keep_cov$V2 == ct,3],split = ',')[[1]]
  sorted_cols <- sort(colSums(obj_it$Y), decreasing = F)
  if(any(!(!(sigs_to_keep_cov_it %in% sorted_cols)))){stop('Some signature name does not match')}
  sigs_to_keep_cov_it <- sigs_to_keep_cov_it[-match(names(sorted_cols[ncol(obj_it$Y)]), sigs_to_keep_cov_it)]
  colnames_all = names(sorted_cols[-ncol(obj_it$Y)])
  diag_cov_names = outer((colnames_all), (colnames_all), function(i,j) paste0(i, '-', j))
  diag_cov_names_lower = diag_cov_names[lower.tri(diag_cov_names)]
  
  colnames_good = as.vector(outer(sigs_to_keep_cov_it, sigs_to_keep_cov_it, function(i,j) paste0(i, '-', j)))
  
  idx_cov_to_infer = match(diag_cov_names_lower, colnames_good)
  idx_cov_to_infer = idx_cov_to_infer[!is.na(idx_cov_to_infer)]
  idx_cov_to_infer
  # return(paste0(ct, "\t", (paste0(idx_cov_to_infer, 'L', collapse=',')), collapse=''))
  return(paste0(ct, "\t", (paste0(idx_cov_to_infer, collapse=',')), collapse=''))
})

write(idx_keep, file = "../../current/subset_sigs_sparse_cov_idx.txt")
