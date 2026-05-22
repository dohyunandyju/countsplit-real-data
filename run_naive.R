run_naive <- function(norm_counts, seed, top.500, pca.scores, null_case = FALSE) {
  set.seed(seed)
  est_cell_type <- kmeans(pca.scores, centers = 2)$cluster # k-means clustering
  
  pvals.naive <- sapply(1:500, function(u) t.test(norm_counts[,top.500][,u]~est_cell_type)$p.value)
  
  pvals.naive_list <- data.frame(Gene = top.500, p_value = pvals.naive)
  # sort p-values in ascending order (lowest at top)
  pvals.naive_list <- pvals.naive_list[order(pvals.naive_list$p_value),]
  
  if (null_case == TRUE) {
    saveRDS(pvals.naive_list, file = paste("pvals.naive.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.naive_list, file = paste("pvals.naive.bcd14.", seed, ".Rds", sep = ""))
  }
  
  print(paste("Naive method for seed ", seed, " completed.", sep = ""))
}