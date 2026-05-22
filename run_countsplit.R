# R script to preprocess data

run_countsplit <- function(raw_counts, seed, null_case = FALSE) {
  
  # Perform countsplit to create train and test sets
  set.seed(seed)
  split_dat <- countsplit(raw_counts)
  Xtrain <- split_dat[[1]]
  Xtest <- split_dat[[2]]
  
  # Pre-process the data (filter for mitochondrial DNA, cells, and features)
  cells.train.seurat <- CreateSeuratObject(counts = Xtrain, min.cells = 3, min.features = 200)
  cells.train.seurat[["percent.mt"]] <- PercentageFeatureSet(cells.train.seurat, pattern = "^MT-")
  cells.train.seurat <- subset(cells.train.seurat, 
                               subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # Modify test set to have same genes as preprocessed training set (so no confusion when testing)
  rows.cells <- rownames(cells.train.seurat)
  cols.cells <- colnames(cells.train.seurat)
  Xtest.subset <- Xtest[rows.cells,cols.cells]
  
  # Normalize and scale data (as well as PCA)
  cells.train.seurat <- NormalizeData(cells.train.seurat)
  cells.train.seurat <- FindVariableFeatures(cells.train.seurat, selection.method = "vst", nfeatures = 2000)
  cells.train.genes <- rownames(cells.train.seurat)
  cells.train.seurat <- ScaleData(cells.train.seurat, features = cells.train.genes)
  cells.train.seurat <- RunPCA(cells.train.seurat, features = VariableFeatures(object = cells.train.seurat))
  
  # Subsetting data (find top 500 most highly variable genes)
  top.500 <- head(VariableFeatures(cells.train.seurat), 500)
  if (null_case == TRUE) {
    saveRDS(top.500, file = paste("top.500.b.cells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(top.500, file = paste("top.500.b.cd14.", seed, ".Rds", sep = ""))
  }
  
  # Clustering
  cells.train.pca_scores <- Embeddings(cells.train.seurat, reduction = "pca")[, 1:10]
  set.seed(seed)
  est_cell_type.countsplit <- kmeans(cells.train.pca_scores, centers=2)$cluster
  
  # Inference (top 500 is top 500 most DE genes found with FindVariableFeatures from Seurat)
  Xtest.t <- t(Xtest.subset / colSums(Xtest.subset)) * 10000 # normalize to size factors and scale
  Xtest.t <- log1p(Xtest.t) # log normalize
  pvals.countsplit <- sapply(1:500, function(u)
    t.test((Xtest.t[,top.500])[,u]~est_cell_type.countsplit)$p.value)
  
  # Need to save this to an RDS object so we can load it into main code
  pvals.countsplit_list <- data.frame(Gene = top.500, p_value = pvals.countsplit)
  # sort p-values in ascending order (lowest at top)
  pvals.countsplit_list <- pvals.countsplit_list[order(pvals.countsplit_list$p_value),]
  
  if (null_case == TRUE) {
    saveRDS(pvals.countsplit_list, file = paste("pvals.countsplit.bcells.", seed, ".Rds", sep = ""))
  } else {
    saveRDS(pvals.countsplit_list, file = paste("pvals.countsplit.bcd14.", seed, ".Rds", sep = ""))
  }
  
  print(paste("Countsplit for seed ", seed, " completed.", sep = ""))
  
}



