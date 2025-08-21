### Here, I provide examples for both Seurat and Harmony integration, briefly discussing the functionality of both approaches.

# Suppose we have read in the data from multiple experiments (biological replicates) of both treatment and control samples, and created objects for each sample
# We have filtered each object separately to remove glial contamination and reserve serotonergic neurons
# We have merged the objects to visualize the relationships between the cells in PC space and UMAP space, but now want to integrate to remove technical variation

# Seurat integration: mutual nearest neighbors (MNNs) and canonical correlation analysis (CCA) to identify anchors or pairs of cells across datasets that capture the same biological state
# Aligns datasets by adjusting expression values of cells such that anchor pairs and their neighbors are matched 

# Starting point: merged_object
# Split merged objectinto a list by batch label or experiment ID
obj.list <- SplitObject(merged_object, split.by = "orig.ident")

# Normalize and find variable features for each dataset
obj.list <- lapply(obj.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)

# Integrate data (creates a new object with "integrated" assay layer)
integrated_object <- IntegrateData(anchorset = anchors, dims = 1:30)

# Standard Seurat workflow on integrated assay
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 50)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 1.0)
integrated <- RunUMAP(integrated, dims = 1:30)

# Visualization
DimPlot(integrated, reduction = "umap", group.by = "orig.ident")

### Harmony integration uses PCA embeddings, not raw expression values like Seurat Integration
# Generates per-batch correction factors and adjusts  PCA coordinates so that cells from different batches align in the same latent space
# Harmony models batch variation statistically as a covariate and removes it in the embedding space

    library(Seurat)
    library(harmony)

# Assuming we are starting with the same unintegrated merged object as above:

# Normalize and find variable features
    merged_object <- NormalizeData(merged_object)
    merged_object <- FindVariableFeatures(merged_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
    merged_object <- ScaleData(merged_object, verbose = FALSE)

# Run PCA
    merged_object <- RunPCA(merged_object, npcs = 50, verbose = FALSE)

# Run Harmony integration: here, "orig.ident" is the batch variable
    merged_object <- RunHarmony(
      object = merged_object,
      group.by.vars = "orig.ident",  
      reduction = "pca",  
      dims.use = 1:50  
    )

# Downstream analysis on Harmony embeddings with harmonized PCA embeddings
    merged_object <- FindNeighbors(merged_object, reduction = "harmony", dims = 1:30)
    merged_object <- FindClusters(merged_object, resolution = 0.5)
    merged_object <- RunUMAP(merged_object, reduction = "harmony", dims = 1:30)
