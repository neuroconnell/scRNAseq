### Here I detail functions that are useful for analyzing data in Seurat

library(dplyr)
library(Seurat)
library(sctransform)
library(pheatmap)
library(edgeR)
library(amap)
library(FactoMineR)
library(scatterD3)
library(ggplots)
library(ggplot2)
library(ClassDiscovery)
library(pheatmap)
library(Matrix)
library(latticeExtra)
library(tools)
library(matrixStats)
library(hdf5r)
library("wesanderson")
library(ape)

# for z-scaled heatmaps
pal = colorRampPalette(c("#1F284C","#2D345D","#6486AD","#D9CCAC","#ECE2C6"))(64) 

# for log2 heatmaps
pal11 = colorRampPalette(c("#2D345D","#6486AD","#D8D8D8","#F9ECE4","#f5d7d0","#f28383","#fc1212"))(64)

# for FeaturePlots
feature_col=colorRampPalette(c("#D0CCAF","#1B288A"))(24)

# colormap for heatmaps in 2020 eLife paper
mycolormap = colorRampPalette(c("royalblue4","ivory","firebrick2"))(256)

# time point colors
time_point_colors_1 =c("#7294D4","#C6CDF7","#E6A0C4","#D8A499")
time_point_colors_2 =c("#7294D4","#C6CDF7","#E6A0C4","#62496F")

orig.ident_colors =colorRampPalette(c("#7294D4","#C6CDF7","#E6A0C4","#D8A499"))(7)

# cluster colors
cluster_colors = colorRampPalette(c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#EBCC2A","#5BBCD6","#046C9A"))(42)
#cluster_colors_2 = colorRampPalette(c("#FF0000","#00A08A","#56bf15","#F2AD00", "#F98400","#EBCC2A","#5BBCD6","#046C9A","#9a0479"))(38)

make_QC_plot <- function(seurat_obj,path){
  fname = paste0(path,deparse(substitute(seurat_obj)),"_QC_plot.pdf")
  seurat_obj[["percent.MT"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribosomal"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS")
  p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.ribosomal"), ncol = 2,pt.size=0)
  ggsave(plot=p,height=8,width=6,dpi=200, filename=fname, useDingbats=FALSE)
}

combine_vlnplots <- function(seurat_list,feature,path){
  fname = paste0(path,feature,"_VlnPlots.pdf")
  # Create an empty list to store the vlnplots
  vlnplots <- list()
  # Loop through the Seurat objects in the list
  for(i in seq_along(seurat_list)){
    vlnplot <- VlnPlot(object = seurat_list[[i]], features = feature,pt.size=0)
    # Add the vlnplot to the list
    vlnplots[[i]] <- vlnplot
  }
  p <- CombinePlots(plots=vlnplots)
  ggsave(plot=p,height=11,width=11,dpi=200, filename=fname, useDingbats=FALSE)
}


save_dendro <- function(seurat_obj,path){
  fname = deparse(substitute(seurat_obj))
  pdf(paste0(path,fname,"_dendro.pdf"),useDingbats=FALSE)
  print(PlotClusterTree(seurat_obj))
  dev.off()
}

save_umap <- function(seurat_obj,group.by=NULL,path,label=TRUE, cols=NULL){	
  fname = deparse(substitute(seurat_obj))
  pdf(paste0(path,fname,"_umap.pdf"),useDingbats=FALSE)
  print(DimPlot(seurat_obj, reduction = "umap",label = label, cols = cols, repel = TRUE,group.by=group.by))
  dev.off()
}

save_FeaturePlot <- function(seurat_obj,features,path){	
  fname = deparse(substitute(seurat_obj))	
  pdf(paste0(path,fname,"_FeaturePlot.pdf"),useDingbats=FALSE)
  print(FeaturePlot(seurat_obj, features=features) &  scale_colour_gradientn(colours = feature_col))
  dev.off()
}

save_VlnPlot <- function(seurat_obj,features,pt.size=0,ncol=1,path){	
  fname = deparse(substitute(seurat_obj))	
  pdf(paste0(path,fname,"_VlnPlot.pdf"),useDingbats=FALSE)
  print(VlnPlot(seurat_obj, features=features,pt.size=pt.size,ncol=ncol))
  dev.off()
}

filter_analysis <- function(seurat_obj,res=1){
  seurat_obj[["percent.MT"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribosomal"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj = subset(seurat_obj, subset = percent.ribosomal > 2 & percent.MT < 15 & nFeature_RNA > 3000 & Fev > 0.5 & Mbp < 2 & Plp1 < 2 & Olig1 == 0 & Aqp4 == 0 & Slc17a6 == 0 & percent.ribosomal > 2)
  return(seurat_obj)
}

standard_analysis <- function(seurat_obj,res=1){
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = c(1:50))
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj = BuildClusterTree(seurat_obj,reorder = TRUE,reorder.numeric = TRUE)
  seurat_obj <- RunUMAP(seurat_obj, dims = c(1:50))
  return(seurat_obj)
}

SCTransform_analysis <- function(seurat_obj,res=1){
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.MT")
  all.genes <- rownames(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = c(1:50))
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  seurat_obj = BuildClusterTree(seurat_obj,reorder = TRUE,reorder.numeric = TRUE)
  seurat_obj <- RunUMAP(seurat_obj, dims = c(1:50))
  return(seurat_obj)
}

return_sizes <- function(seurat_list){
  for(i in seq_along(seurat_list)){
    print(dim(seurat_list[[i]]))
  }
}

create_folder <- function(folder_name) {
  # Check if the folder name is valid
  if (!is.character(folder_name) || length(folder_name) != 1) {
    stop("Invalid folder name")
  }
  # Check if the folder already exists
  if (dir.exists(folder_name)) {
    message("Folder already exists")
  } else {
    # Create the folder using dir.create function
    dir.create(folder_name)
    message("Folder created successfully")
  }
}

find_markers_heatmap <- function(seurat_obj, only.pos=TRUE,path) {
  create_folder(path)
  # Find all markers of cluster identity
  markers <- FindAllMarkers(seurat_obj, only.pos = only.pos, min.pct = 0.25, logfc.threshold = 0.25)
  # Group markers by cluster
  markers_by_cluster <- split(markers, markers$cluster)
  # Loop over each cluster
  for (i in 1:length(markers_by_cluster)) {
    # Get the top 50 markers for the current cluster
    top_markers <- head(markers_by_cluster[[i]], n = 50)$gene
    # Create a heatmap of the top markers using DoHeatmap
    pdf(paste0(path,"cluster_",i,"_heatmap.pdf"),useDingbats=FALSE)
    print(DoHeatmap(seurat_obj, features = top_markers, group.by = "ident", label = FALSE, assay="RNA") + scale_fill_gradientn(colors = pal))
    dev.off()
  }
}

output_markers_heatmap_integrated <- function(seurat_obj, markers, only.pos=TRUE,path,group.colors=NULL) {
  create_folder(path)
  # Group markers by cluster
  markers_by_cluster <- split(markers, markers$cluster)
  # Loop over each cluster
  for (i in 1:length(markers_by_cluster)) {
    # Get the top 50 markers for the current cluster
    top_markers <- head(markers_by_cluster[[i]], n = 50)$gene
    # Create a heatmap of the top markers using DoHeatmap
    pdf(paste0(path,"cluster_",i,"_heatmap.pdf"),useDingbats=FALSE)
    print(DoHeatmap(seurat_obj, features = top_markers, group.by = "ident", group.colors = group.colors, label = FALSE) + scale_fill_gradientn(colors = pal))
    dev.off()
  }
}


output_markers_heatmap_RNA <- function(seurat_obj, markers, path,group.colors=NULL) {
  create_folder(path)
  # Group markers by cluster
  markers_by_cluster <- split(markers, markers$cluster)
  # Loop over each cluster
  for (i in 1:length(markers_by_cluster)) {
    # Get the top 50 markers for the current cluster
    top_markers <- head(markers_by_cluster[[i]], n = 50)$gene
    # Create a heatmap of the top markers using DoHeatmap
    seurat_obj <- ScaleData(seurat_obj, features = top_markers)
    pdf(paste0(path,"cluster_",i,"_heatmap.pdf"),useDingbats=FALSE)
    print(DoHeatmap(seurat_obj, features = top_markers, group.by = "ident", group.colors = group.colors, label = FALSE) + scale_fill_gradientn(colors = pal))
    dev.off()
  }
}

# Define the function
stacked_barplot <- function(seurat_obj, meta_var1, meta_var2, colors = NULL, orientation = "vertical", reverse_y_order = FALSE) {
  # Check if metadata variables exist in the Seurat object
  if (!(meta_var1 %in% colnames(seurat_obj@meta.data)) | !(meta_var2 %in% colnames(seurat_obj@meta.data))) {
    stop("One or both metadata variables do not exist in the Seurat object.")
  }
  
  # Validate the orientation argument
  if (!orientation %in% c("vertical", "horizontal")) {
    stop("Invalid orientation value. Please use 'vertical' or 'horizontal'.")
  }
  
  # Extract the metadata
  metadata <- seurat_obj@meta.data
  
  # Convert both meta_var1 and meta_var2 to factors to ensure proper handling in ggplot
  metadata[[meta_var1]] <- as.factor(metadata[[meta_var1]])
  metadata[[meta_var2]] <- as.factor(metadata[[meta_var2]])
  
  # Reverse the order of meta_var2 if reverse_y_order is TRUE and orientation is horizontal
  if (orientation == "horizontal" && reverse_y_order) {
    metadata[[meta_var2]] <- factor(metadata[[meta_var2]], levels = rev(levels(metadata[[meta_var2]])))
  }
  
  # Create a contingency table and convert it to long format for ggplot2
  prop_table <- metadata %>%
    group_by(!!sym(meta_var2), !!sym(meta_var1)) %>%
    summarise(n = n()) %>%
    mutate(percentage = n / sum(n) * 100)
  
  # Generate dynamic colors if not provided
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(unique(metadata[[meta_var1]])))
  } else {
    if (length(colors) != length(unique(metadata[[meta_var1]]))) {
      stop("Number of colors provided does not match the number of levels in meta_var1.")
    }
  }
  
  # Create the ggplot object
  p <- ggplot(prop_table, aes_string(x = meta_var2, y = "percentage", fill = meta_var1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    ggtitle(paste("Percent composition of", meta_var1, "in terms of", meta_var2))
  
  # Adjust the plot based on the orientation
  if (orientation == "vertical") {
    p <- p + 
      xlab(meta_var2) +
      ylab("Percentage") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Angle x-axis labels
      scale_x_discrete(limits = levels(metadata[[meta_var2]]))  # Explicitly set the levels
  } else if (orientation == "horizontal") {
    p <- p + 
      xlab("Percentage") +  # Label x-axis for percentage
      ylab(meta_var2) +  # Use meta_var2 as the label for the y-axis
      coord_flip() +  # Flip coordinates for horizontal orientation
      scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Format y-axis as percentage after flipping
      scale_x_discrete(limits = levels(metadata[[meta_var2]])) +  # Ensure y-axis (original x-axis) is treated as discrete
      theme(axis.text.y = element_text(angle = 45, hjust = 1))  # Angle y-axis labels
  }
  
  return(p)
}

PlotLogNormalizedHistogram <- function(seurat_obj, bins = 50, log_expression_range = NULL, plot_mean = FALSE, highlight_gene = NULL, cell_identity = NULL) {
  # Check if the input is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  # Subset Seurat object based on cell identities if specified
  if (!is.null(cell_identity)) {
    cell_identity <- as.character(cell_identity)
    missing_idents <- setdiff(cell_identity, levels(Idents(seurat_obj)))
    if (length(missing_idents) > 0) {
      stop(paste("The following identities are not found in the Seurat object:", paste(missing_idents, collapse = ", ")))
    }
    seurat_obj <- subset(seurat_obj, idents = cell_identity)
  }
  
  # Extract log-normalized expression data
  log_norm_data <- GetAssayData(seurat_obj, slot = "data", assay = "RNA")
  
  # Calculate mean values across cells if requested
  if (plot_mean) {
    log_norm_values <- rowMeans(log_norm_data)
  } else {
    log_norm_values <- as.vector(log_norm_data)
  }
  
  # Filter values if a range is specified
  if (!is.null(log_expression_range)) {
    log_norm_values <- log_norm_values[log_norm_values >= log_expression_range[1] &
                                         log_norm_values <= log_expression_range[2]]
  }
  
  # Ensure log_norm_values is converted into a data frame for ggplot
  hist_data <- data.frame(log_expression = log_norm_values)
  
  # Prepare gene highlights
  highlight_positions <- NULL
  highlight_labels <- NULL
  if (!is.null(highlight_gene)) {
    # Ensure highlight_gene is a character vector
    highlight_gene <- as.character(highlight_gene)
    
    # Check that all genes exist in the Seurat object
    missing_genes <- setdiff(highlight_gene, rownames(log_norm_data))
    if (length(missing_genes) > 0) {
      stop(paste("The following genes were not found in the Seurat object:", paste(missing_genes, collapse = ", ")))
    }
    
    # Calculate positions for the highlighted genes
    if (plot_mean) {
      highlight_positions <- sapply(highlight_gene, function(g) mean(log_norm_data[g, ]))
    } else {
      highlight_positions <- sapply(highlight_gene, function(g) as.vector(log_norm_data[g, ]))
    }
    highlight_labels <- highlight_gene
  }
  
  # Create the histogram
  title_identity <- ifelse(is.null(cell_identity), "All Cells", paste(cell_identity, collapse = ", "))
  p <- ggplot(hist_data, aes(x = log_expression)) +
    geom_histogram(bins = bins, fill = "blue", alpha = 0.7) +
    theme_minimal() +
    labs(title = if (plot_mean) paste("Histogram of Mean Log-Normalized Expression Levels for", title_identity) else paste("Histogram of Log-Normalized Expression Levels for", title_identity),
         x = if (plot_mean) "Mean Log-Normalized Expression" else "Log-Normalized Expression",
         y = "Frequency") +
    theme(
      text = element_text(size = 14)
    )
  
  # Add vertical lines and labels for the highlighted genes
  if (!is.null(highlight_gene) && !is.null(highlight_positions)) {
    highlight_data <- data.frame(position = highlight_positions, label = highlight_labels)
    
    p <- p +
      geom_vline(data = highlight_data, aes(xintercept = position), color = "red", linetype = "dashed", size = 1) +
      geom_text(data = highlight_data, aes(x = position, 
                                           y = max(table(cut(log_norm_values, bins))) * 0.6, # Position slightly below
                                           label = label), 
                angle = 45, color = "red", hjust = -0.1, vjust = 1, size = 4)
  }
  
  return(p)
}

CalculateGeneSetScore <- function(seurat_obj, gene_list, score_name = "GeneSetScore") {
  # Ensure all genes in the list are in the dataset
  valid_genes <- intersect(gene_list, rownames(seurat_obj))
  if (length(valid_genes) == 0) {
    stop("None of the provided genes are found in the Seurat object.")
  }
  
  # Get the data matrix for valid genes
  gene_matrix <- GetAssayData(seurat_obj, slot = "data")[valid_genes, , drop = FALSE]
  
  # Normalize expression values to range 0-1 for each gene
  normalized_gene_matrix <- apply(gene_matrix, 1, function(x) {
    if (max(x) == min(x)) {
      # Handle cases where all values are the same (to avoid division by zero)
      return(rep(0, length(x)))
    } else {
      return((x - min(x)) / (max(x) - min(x)))
    }
  })
  
  # Ensure output is a matrix (even if only one gene is present)
  if (is.vector(normalized_gene_matrix)) {
    normalized_gene_matrix <- t(normalized_gene_matrix)
  }
  
  # Sum normalized expression values across genes for each cell
  cell_scores <- rowSums(normalized_gene_matrix)
  
  # Add the scores to the Seurat object's metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = cell_scores, col.name = score_name)
  
  return(seurat_obj)
}

CalculateGeneRatio <- function(seurat_obj, gene1, gene2, ratio_name = "GeneRatio") {
  # Ensure both genes are present in the Seurat object
  valid_genes <- c(gene1, gene2) %in% rownames(seurat_obj)
  if (!all(valid_genes)) {
    stop(paste("One or both genes are not found in the Seurat object:",
               paste(c(gene1, gene2)[!valid_genes], collapse = ", ")))
  }
  
  # Extract normalized expression values for the genes
  gene1_expr <- GetAssayData(seurat_obj, slot = "data")[gene1, ]
  gene2_expr <- GetAssayData(seurat_obj, slot = "data")[gene2, ]
  
  # Calculate the ratio (gene1 / gene2)
  gene_ratio <- gene1_expr / gene2_expr
  
  # Handle cases where gene2_expr is zero (to avoid division by zero)
  gene_ratio[is.infinite(gene_ratio) | is.nan(gene_ratio)] <- NA
  
  # Add the ratio to the Seurat object's metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = gene_ratio, col.name = ratio_name)
  
  return(seurat_obj)
}

CountCellsByGroup <- function(seurat_obj, group_label) {
  # Ensure the grouping label is valid
  if (!group_label %in% colnames(seurat_obj@meta.data)) {
    stop(paste("The group label", group_label, "is not found in the metadata of the Seurat object."))
  }
  
  # Extract the grouping information
  group_info <- seurat_obj@meta.data[[group_label]]
  
  # Count cells in each group
  cell_counts <- table(group_info)
  
  # Return the counts as a named vector
  return(as.data.frame(cell_counts, stringsAsFactors = FALSE))
}



CalculateGeneSetScore <- function(seurat_obj, gene_list, score_name = "GeneSetScore") {
  # Ensure all genes in the list are in the dataset
  valid_genes <- intersect(gene_list, rownames(seurat_obj))
  if (length(valid_genes) == 0) {
    stop("None of the provided genes are found in the Seurat object.")
  }
  
  # Get the data matrix for valid genes
  gene_matrix <- GetAssayData(seurat_obj, slot = "data")[valid_genes, , drop = FALSE]
  
  # Normalize expression values to range 0-1 for each gene
  normalized_gene_matrix <- apply(gene_matrix, 1, function(x) {
    if (max(x) == min(x)) {
      # Handle cases where all values are the same (to avoid division by zero)
      return(rep(0, length(x)))
    } else {
      return((x - min(x)) / (max(x) - min(x)))
    }
  })
  
  # Ensure output is a matrix (even if only one gene is present)
  if (is.vector(normalized_gene_matrix)) {
    normalized_gene_matrix <- t(normalized_gene_matrix)
  }
  
  # Sum normalized expression values across genes for each cell
  cell_scores <- rowSums(normalized_gene_matrix)
  
  # Add the scores to the Seurat object's metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = cell_scores, col.name = score_name)
  
  return(seurat_obj)
}

calculate_grouping_skewness <- function(seurat_obj, group_var1, group_var2, method = c("variance", "cv", "gini", "entropy"), plot = FALSE) {
  # Extract metadata
  meta_data <- seurat_obj@meta.data
  
  # Ensure grouping variables exist in the metadata
  if (!(group_var1 %in% colnames(meta_data)) || !(group_var2 %in% colnames(meta_data))) {
    stop("Specified grouping variables must exist in the Seurat object's metadata.")
  }
  
  # Ensure valid method
  method <- match.arg(method)
  
  # Calculate counts for each combination of group_var1 and group_var2
  counts_table <- table(meta_data[[group_var1]], meta_data[[group_var2]])
  
  # Convert to proportions within each level of group_var1
  prop_table <- prop.table(counts_table, margin = 1)
  
  # Function to calculate Shannon entropy
  calculate_entropy <- function(proportions) {
    proportions <- proportions[proportions > 0] # Remove zeroes to avoid log issues
    -sum(proportions * log2(proportions))
  }
  
  # Function to calculate coefficient of variation
  calculate_cv <- function(proportions) {
    sd(proportions) / mean(proportions)
  }
  
  # Function to calculate Gini coefficient
  calculate_gini <- function(proportions) {
    gini(proportions)
  }
  
  # Apply the chosen method
  skewness_per_group <- switch(method,
                               "variance" = apply(prop_table, 1, var),
                               "cv" = apply(prop_table, 1, calculate_cv),
                               "gini" = apply(prop_table, 1, calculate_gini),
                               "entropy" = apply(prop_table, 1, calculate_entropy)
  )
  
  # Format the result as a data frame for easy interpretation
  result <- data.frame(
    group_var1_level = rownames(prop_table),
    skewness = skewness_per_group
  )
  
  # Convert group_var1_level to numeric if possible and reorder
  result$group_var1_level <- as.numeric(as.character(result$group_var1_level))
  result <- result[order(result$group_var1_level), ]
  
  # Optionally plot the result
  if (plot) {
    plot <- ggplot(result, aes(x = factor(group_var1_level), y = skewness)) +
      geom_bar(stat = "identity", fill = "steelblue", color = "black") +
      scale_x_discrete(labels = result$group_var1_level) +
      labs(
        title = paste("Skewness of", group_var2, "in", group_var1, "Clusters"),
        x = group_var1,
        y = paste("Skewness (", method, ")", sep = "")
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_line(color = "gray", linetype = "dotted")
      )
    print(plot)
  }
  
  return(result)
}

exclude_number <- function(exclude, max_number) {
  return(setdiff(1:max_number, exclude))
}
