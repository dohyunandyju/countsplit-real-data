run_cadet_null <- function(norm_counts, seed, top.500, null_case = FALSE) {
  
  # Calculate the covariance matrix (without demean)
  cov_null <- coop::covar(norm_counts[,top.500])
  
  # CADET null case with non-zero off diagonals
  pvals.cadet.null <- sapply(1:500, function(u) 
    kmeans_inference_1f(norm_counts[,top.500], k = 2, 1, 2, feat = u, 
                        iso = FALSE, sig = NULL, covMat = cov_null, seed = seed)$pval)
  
  pvals.cadet.null_list <- data.frame(Gene = top.500, p_value = pvals.cadet.null)
  # sort p-values in ascending order (lowest at top)
  pvals.cadet.null_list <- pvals.cadet.null_list[order(pvals.cadet.null_list$p_value),]
  
  if (null_case == TRUE) {
    saveRDS(pvals.cadet.null_list, file = paste("pvals.cadet.null.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.cadet.null_list, file = paste("pvals.cadet.null.bcd14.", seed, ".Rds", sep = ""))
  }
  
  # CADET null case with zero off diagonals
  cov_indep <- cov_null
  cov_indep[row(cov_indep) != col(cov_indep)] <- 0
  
  pvals.cadet.indep_null <- sapply(1:500, function(u) 
    kmeans_inference_1f(norm_counts[,top.500], k = 2, 1, 2, feat = u, 
                        iso = FALSE, sig = NULL, covMat = cov_indep, seed = seed)$pval)
  
  pvals.cadet.indep_null_list <- data.frame(Gene = top.500, p_value = pvals.cadet.indep_null)
  # sort p-values in ascending order (lowest at top)
  pvals.cadet.indep_null_list <- pvals.cadet.indep_null_list[order(pvals.cadet.indep_null_list$p_value),]

  if (null_case == TRUE) {
    saveRDS(pvals.cadet.indep_null_list, file = paste("pvals.cadet.indep_null.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.cadet.indep_null_list, file = paste("pvals.cadet.indep_null.bcd14.", seed, ".Rds", sep = ""))
  }
  
  print(paste("CADET null for seed ", seed, " completed.", sep = ""))
}