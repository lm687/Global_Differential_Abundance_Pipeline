## for each cancer type under consideration, write


sigs_to_keep_cov <- read.table("../../current/subset_sigs_sparse_cov.txt", sep = "\t")
sigs_to_keep_cov <- sigs_to_keep_cov[!(sigs_to_keep_cov$V3 == "?"),]

sigs_to_keep_cov_nonexo <- read.table("../../current/subset_sigs_sparse_cov_nonexo.txt", sep = "\t")
sigs_to_keep_cov_nonexo <- sigs_to_keep_cov_nonexo[!(sigs_to_keep_cov_nonexo$V3 == "?"),]

give_indices <- function(ct, remove_exogenous=F, sigs_to_keep_cov_df){
  if(!(remove_exogenous)){
    obj_it = load_PCAWG(ct, typedata = "signatures")
  }else{
    nonexogenous = read.table("../../data/cosmic/exogenous_signatures_SBS.txt", sep = "\t", comment.char = "#", fill = F)
    obj_it <- give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = "signatures"),
                                          sigs_to_remove = unique(nonexogenous$V1))
  }
  sigs_to_keep_cov_it <- strsplit(sigs_to_keep_cov_df[sigs_to_keep_cov_df$V2 == ct,3],split = ',')[[1]]
  sorted_cols <- sort(colSums(obj_it$Y), decreasing = F)
  if(any(!(!(sigs_to_keep_cov_it %in% sorted_cols)))){stop('Some signature name does not match')}
  sigs_to_keep_cov_it <- sigs_to_keep_cov_it[-match(names(sorted_cols[ncol(obj_it$Y)]), sigs_to_keep_cov_it)]
  colnames_all = names(sorted_cols[-ncol(obj_it$Y)])
  diag_cov_names = (outer((colnames_all), (colnames_all), function(i,j) paste0(i, '-', j)))
  diag_cov_names_lower = diag_cov_names[lower.tri(diag_cov_names)]
  
  colnames_good = as.vector(outer(sigs_to_keep_cov_it, sigs_to_keep_cov_it, function(i,j) paste0(i, '-', j)))
  
  idx_cov_to_infer = match(colnames_good, diag_cov_names_lower)
  # .mat_reindexing = matrix(NA, ncol=(ncol(obj_it$Y)-1), nrow=(ncol(obj_it$Y)-1))
  # matrix(1:((ncol(obj_it$Y)-1)**2), ncol=(ncol(obj_it$Y)-1), nrow=(ncol(obj_it$Y)-1))
  # .mat_reindexing[lower.tri(.mat_reindexing)] = 1:( ((ncol(obj_it$Y)-1)**2 - (ncol(obj_it$Y)-1)) /2)
  # idx_cov_to_infer =.mat_reindexing[idx_cov_to_infer]
  
  idx_cov_to_infer = idx_cov_to_infer[!is.na(idx_cov_to_infer)]
  idx_cov_to_infer
  # return(paste0(ct, "\t", (paste0(idx_cov_to_infer, 'L', collapse=',')), collapse=''))
  return(paste0(ct, "\t", (paste0(idx_cov_to_infer, collapse=',')), collapse=''))
}

idx_keep = sapply(sigs_to_keep_cov$V2, give_indices, sigs_to_keep_cov_df=sigs_to_keep_cov)
write("## These are the indices in the covariance matrix which we are going to keep", file = "../../current/subset_sigs_sparse_cov_idx.txt",append = F)
write(idx_keep, file = "../../current/subset_sigs_sparse_cov_idx.txt",append = T)


# idx_keep_nonexo = sapply(sigs_to_keep_cov_nonexo$V2, give_indices, remove_exogenous=T, sigs_to_keep_cov_df=sigs_to_keep_cov_nonexo)
idx_keep_nonexo = sapply(sigs_to_keep_cov_nonexo$V2, give_indices, remove_exogenous=T, sigs_to_keep_cov_df=sigs_to_keep_cov_nonexo)
write("## These are the indices in the covariance matrix which we are going to keep", file = "../../current/subset_sigs_sparse_cov_idx_nonexo.txt",append = F)
write(idx_keep_nonexo, file = "../../current/subset_sigs_sparse_cov_idx_nonexo.txt",append = T)
