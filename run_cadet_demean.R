run_cadet_demean <- function(norm_counts, seed, top.500, null_case = FALSE) {
  # create the demean matrix
  set.seed(seed)
  kmeans_clust <- kmeans_estimation(norm_counts[,top.500], k = 2, iter.max = 30, seed = seed)
  final_clust <- kmeans_clust$cluster[[kmeans_clust$iter]]
  final_centroids <- kmeans_clust$centers[[kmeans_clust$iter]]
  counts_demean <- norm_counts[,top.500] - final_centroids[final_clust, , drop = FALSE]
  cov_demean <- coop::covar(counts_demean)
  
  # CADET demean case with non-zero off diagonals
  pvals.cadet.demean <- sapply(1:500, function(u) 
    kmeans_inference_1f(norm_counts[,top.500], k = 2, 1, 2, feat = u, 
                        iso = FALSE, sig = NULL, covMat = cov_demean, seed = seed)$pval)
  
  pvals.cadet.demean_list <- data.frame(Gene = top.500, p_value = pvals.cadet.demean)
  # sort p-values in ascending order (lowest at top)
  pvals.cadet.demean_list <- pvals.cadet.demean_list[order(pvals.cadet.demean_list$p_value),]
  
  if (null_case == TRUE) {
    saveRDS(pvals.cadet.demean_list, file = paste("pvals.cadet.demean.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.cadet.demean_list, file = paste("pvals.cadet.demean.bcd14.", seed, ".Rds", sep = ""))
  }
  
  # CADET demean case with zero off diagonals
  cov_indep_demean <- cov_demean
  cov_indep_demean[row(cov_indep_demean) != col(cov_indep_demean)] <- 0
  
  pvals.cadet.indep_demean <- sapply(1:500, function(u) 
    kmeans_inference_1f(norm_counts[,top.500], k = 2, 1, 2, feat = u, 
                        iso = FALSE, sig = NULL, covMat = cov_indep_demean, seed = seed)$pval)
  
  pvals.cadet.indep_demean_list <- data.frame(Gene = top.500, p_value = pvals.cadet.indep_demean)
  # sort p-values in ascending order (lowest at top)
  pvals.cadet.indep_demean_list <- pvals.cadet.indep_demean_list[order(pvals.cadet.indep_demean_list$p_value),]
  
  if (null_case == TRUE) {
    saveRDS(pvals.cadet.indep_demean_list, file = paste("pvals.cadet.indep_demean.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.cadet.indep_demean_list, file = paste("pvals.cadet.indep_demean.bcd14.", seed, ".Rds", sep = ""))
  }
  
  print(paste("CADET demean for seed ", seed, " completed.", sep = ""))
}