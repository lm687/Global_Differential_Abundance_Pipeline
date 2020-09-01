duplicate_rows = function(df, n){
  .x <- do.call('rbind', lapply(1:nrow(df), function(i) (do.call('rbind', replicate(n, df[i,], simplify=FALSE )))))
  #rownames(x) <- paste0(1:(nrow(.x)))
  .x
}

# softmax = function(v){
#   .x = sapply(v, exp)
#   .x/sum(.x)
# }
# 
# softmax_mat = function(v){
#   .x = t(apply(v, 1, exp))
#   sweep(.x, 1, rowSums(.x), '/')
# }

normalise_rw = function (x) 
{
  sweep(x, 1, rowSums(x), "/")
}

