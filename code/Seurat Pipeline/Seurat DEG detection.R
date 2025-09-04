### Here, we will be predicting differentially expressed genes. We will be testing for differential expression of genes due to treatment or cell-identity. 
# The native Seurat function for predicting differentially expressed genes between clusters is "FindAllMarkers".
# The native Seurat function for predicting differentially expressed genes between identity classes (e.g., repliclates) is "FindMarkers".

### FindAllMarkers: Identfying differentially expressed genes between cell populations
# Here, we assume that we have a Seurat object containing purified serotonergic neurons, with 2 classes, treatment and control

# 1. Set orig.ident as the active identity --> level of resolution for this analysis are cell clusters 
Idents(serotonin_obj) <- "seurat_clusters"

# 2. If the data has been scaled, integrated, or parsed into separate layers, layers must be joined
serotonin_obj_joined<-JoinLayers(Serotonin_Obj)

# Confirm the identity levels (these represent your samples/runs)
levels(serotonin_obj_joined)

#Run FindAllMarkers to predict differentially expressed genes between cell populations (clusters)
#Example: Compare one sample (e.g., Sample1) vs. all others
serotonin_cluster_markers <- FindAllMarkers(
  object = serotonin_obj_joined,
  only.pos = TRUE,            # Only return positive markers
  min.pct = 0.25,             # Minimum % of cells expressing gene in either group
  logfc.threshold = 0.5,      # Minimum log2 fold change
  assay = "RNA"               # Specify assay if needed
)

# Write the data to an excel file for analysis
write.csv(serotonin_cluster_markers, file = "Serotonergic_Subpopulation_Markers", row.names = T)


### FindMarkers: Identifying differentially expressed genes across classes (e.g., replicates)
# We assume again that we are starting with a Seurat Object containing purified serotonergic neurons with 2 classes, treatment and control

# 1. Here, we set the active identity to the level of the class, which is stored as "orig.ident", 
Idents(serotonin_obj) <- "orig.ident"

# 2. If the data has been scaled, integrated, or parsed into separate layers, layers must be joined
serotonin_obj_joined<-JoinLayers(serotonin_obj)

# Confirm the identity levels (these represent your samples/runs)
levels(serotonin_obj_joined)

#Run FindMarkers to predict differentially expressed genes between classes (replicates)
#Example: Compare one sample (e.g., Sample1) vs. all others
serotonin_cluster_markers <- FindMarkers(
  object = serotonin_obj_joined,
  only.pos = TRUE,            # Only return positive markers
  min.pct = 0.25,             # Minimum % of cells expressing gene in either group
  logfc.threshold = 0.5,      # Minimum log2 fold change
  assay = "RNA"               # Specify assay if needed
)

# Write the data to an excel file for analysis
write.csv(serotonin_cluster_markers, file = "Replicate_Gene_Markers", row.names = T)




### We likely want to identify top DEGs at the level of condition, not between all subjects. For this, we will create a new layer called "condition". 
#In this layer, we will assign cells from the respective replicate (orig.ident) and then test for DEGs using FindMarkers

# Define treatment and control replicates
treatment_reps <- c("treated_replicate_1", "treated_replicate_2")
control_reps   <- c("control_replicate_1", "control_replicate_2")

# Add a 'condition' column to metadata
serotonin_obj$condition <- case_when(
  serotonin_obj$orig.ident %in% treatment_reps ~ "Treatment",
  serotonin_obj$orig.ident %in% control_reps   ~ "Control"
)

# Convert to factor
serotonin_obj$condition <- factor(serotonin_obj$condition, levels = c("Control", "Treatment"))

# If the data has been scaled, integrated, or parsed into separate layers, layers must be joined
serotonin_obj_joined<-JoinLayers(serotonin_obj)
# Set identities to 'condition'
Idents(serotonin_obj_joined) <- "condition"

# Run FindMarkers at condition level
condition_markers <- FindMarkers(
  object       = serotonin_obj_joined,
  ident.1      = "Treatment",
  ident.2      = "Control",
  only.pos     = TRUE,
  min.pct      = 0.25,
  logfc.threshold = 0.5,
  assay        = "RNA"
)

# Save to CSV
write.csv(condition_markers, file = "Condition_DEGs.csv", row.names = TRUE)























