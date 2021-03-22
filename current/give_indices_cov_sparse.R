## for each cancer type under consideration, write


sigs_to_keep_cov = read.table("../../current/subset_sigs_sparse_cov.txt", sep = "\t")
sigs_to_keep_cov <- strsplit(sigs_to_keep_cov[sigs_to_keep_cov$V1 == "CNS_Medullo",2],split = ',')[[1]]
sorted_cols <- sort(colSums(CNS_Medullo$Y), decreasing = F)
if(any(!(!(sigs_to_keep_cov %in% sorted_cols)))){stop('Some signature name does not match')}
sigs_to_keep_cov <- sigs_to_keep_cov[-match(names(sorted_cols[ncol(CNS_Medullo$Y)]), sigs_to_keep_cov)]
colnames_all = names(sorted_cols[-ncol(CNS_Medullo$Y)])
diag_cov_names = outer((colnames_all), (colnames_all), function(i,j) paste0(i, '-', j))
diag_cov_names_lower = diag_cov_names[lower.tri(diag_cov_names)]

colnames_good = as.vector(outer(sigs_to_keep_cov, sigs_to_keep_cov, function(i,j) paste0(i, '-', j)))

idx_cov_to_infer = match(diag_cov_names_lower, colnames_good)
idx_cov_to_infer = idx_cov_to_infer[!is.na(idx_cov_to_infer)]
idx_cov_to_infer
cat(paste0(idx_cov_to_infer, 'L', sep=','))

