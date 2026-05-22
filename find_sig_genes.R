find_sig_genes <- function(pvals_list, p_val) {
  newlist <- c()
  
  for(list in pvals_list) {
    sub_list <- subset(list, p_value < p_val)
    newlist <- c(newlist, nrow(sub_list))
  }
  return(newlist)
}