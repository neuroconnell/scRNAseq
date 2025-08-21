### Read in .h5 matrix
Example = Read10X_h5("Example_filtered_feature_bc_matrix.h5")

### Create Seurat object and conduct QC for mitochondrial and ribosomal genes
Example_Seurat_Object <- CreateSeuratObject(counts = Example, project = "Example", min.cells = 3, min.features = 500)
Example_Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Example_Seurat_Object, pattern = "^mt-")
Example_Seurat_Object[["percent.ribosomal"]] <- PercentageFeatureSet(Example_Seurat_Object, pattern = "^Rps")

### Visualize QC metrics
VlnPlot(Example_Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribosomal"), ncol = 4, pt.size=0)
## Subset object by QC cutoffs
Example_Seurat_Object <- subset(Example_Seurat_Object, subset = percent.mt < 15 & nFeature_RNA < 10000 & nFeature_RNA > 500)

### Normalize, find variable features, and scale data
Example_Seurat_Object <- NormalizeData(Example_Seurat_Object, normalization.method = "LogNormalize", scale.factor = 10000)
Example_Seurat_Object <- FindVariableFeatures(Example_Seurat_Object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Example_Seurat_Object)
Example_Seurat_Object <- ScaleData(Example_Seurat_Object, features = all.genes)

### Run principal component analysis, find neighbors, and Louvain clustering
Example_Seurat_Object <- RunPCA(Example_Seurat_Object, features = VariableFeatures(object = Example_Seurat_Object))
Example_Seurat_Object <- FindNeighbors(Example_Seurat_Object, dims = c(1:50))
Example_Seurat_Object <- FindClusters(Example_Seurat_Object, resolution = 1.0)

### Generate UMAP plot using first 50 PCs
Example_Seurat_Object <- RunUMAP(Example_Seurat_Object, dims = c(1:50))

### Define cluster colors
cluster_colors = colorRampPalette(c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#EBCC2A","#5BBCD6","#046C9A"))(39)

### Print DimPlot
print(DimPlot(Example_Seurat_Object, reduction = "umap", label = TRUE, repel = TRUE, cols=cluster_colors))

### Evaluate nFeature_RNA by cluster
FeaturePlot(Example_Seurat_Object, features = c("nFeature_RNA"))
VlnPlot(Example_Seurat_Object, features = "nFeature_RNA")

### Visualize violin plot or feature plot of cluster-specific expression for known serotonergic markers
VlnPlot(Example_Seurat_Object, features=c("Fev","Lmx1b","Tph2","Slc6a4"), pt.size=0.1, cols=cluster_colors)
FeaturePlot(Example_Seurat_Object, features= c("Fev","Lmx1b","Tph2","Slc6a4"))


