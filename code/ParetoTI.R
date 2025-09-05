### Here, I use Pareto task inference (ParetoTI) to identify archetypal gene expression states by bounding the cells embedded in the n-dimensional principal component space with an m-dimensional convex polytope. 

# ParetoTI seeks to understand the dynamics of transcriptomic "prioritization" of a given cell population, that is, how specialized a particular cell is.
# The vertices of the m-dimensional shape comprise "archetypal" gene expression states; the location of every cell within the volume of the shape are weighted combinations (alphas) of these vertices. 
# Cells that occupy faces, edges, or vertices are considered "specialists" whereas cells that occupy the volume may be considered "generalists". 

# I calculate the number of vertices (and edges/faces) by iteratively optimizing specific parameters:
  # 1. Appropriately fit the correct number of vertices to capture the embedded cells (residual sum of squares)
  # 2. Assign archetypessuch that each archetypes comprise interpretable partions (silhouette score)

# I then visualize the contributions (weights/alphas) of each vertex to the gene expression state of a given cell, visualized using violin plots and feature plots.

# Next, I depict the m-dimensional polytope in 3 dimensions with an interactable browser, and identify the gene expression states associated with each archetype.

###ParetoTI Workflo
# Required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(archetypes)
library(gprofiler2)
library(cluster)
library(plotly)
library(scales)

# Example cluster colors
cluster_colors = colorRampPalette(
  c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#EBCC2A","#5BBCD6","#046C9A")
)(30)

# DimPlot for visualization
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = cluster_colors) +
  ggtitle('All cells') +
  theme(plot.title = element_text(hjust = 0.5))

# Extract PCA coordinates
DefaultAssay(seurat_obj) <- "RNA"
pca_coords <- Embeddings(seurat_obj, "pca")[, 1:5]

# Determine optimal k (elbow plot with RSS)
rss_values <- sapply(3:15, function(k) {
  model <- stepArchetypes(pca_coords, k = k, nrep = 10, verbose = FALSE)
  best <- bestModel(model)
  best$rss
})

plot(3:15, rss_values, type = 'b', pch = 19,
     xlab = "Number of Archetypes (k)", ylab = "Residual Sum of Squares",
     main = "Elbow Plot for Archetype Selection")

# Final model with chosen k
k <- 12
final_model <- bestModel(stepArchetypes(pca_coords, k = k, nrep = 5, verbose = FALSE))

# Silhouette scores (cluster separation quality)
alphas <- final_model$alphas
dominant_archetype <- apply(alphas, 1, which.max)
dominant_archetype <- factor(dominant_archetype, labels = paste0("A", 1:ncol(alphas)))

# Reduce alpha space for distance computation
pca_res <- prcomp(alphas, center = TRUE, scale. = FALSE)
coords <- pca_res$x[, 1:min(10, ncol(pca_res$x))]

# Compute silhouette scores
dists <- dist(coords)
sil <- silhouette(as.integer(dominant_archetype), dists)
sil_df <- as.data.frame(sil)
sil_df$archetype <- dominant_archetype[as.numeric(rownames(sil_df))]

# Plot silhouette scores
ggplot(sil_df, aes(x = archetype, y = sil_width, fill = archetype)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2, size = 0.5) +
  theme_minimal() +
  labs(title = "Silhouette Scores by Dominant Archetype",
       x = "Archetype", y = "Silhouette Width") +
  theme(legend.position = "none")

# Interactive 3D visualization with cluster & archetype sliders
cluster_labels <- factor(seurat_obj@meta.data$cluster_id)  # replace with your cluster column
proj <- final_model$alphas
pca_res <- prcomp(proj, center = TRUE, scale. = FALSE)
pca_coords <- pca_res$x[, 1:3]
colnames(pca_coords) <- c("PC1", "PC2", "PC3")

cells_df <- data.frame(
  PC1 = pca_coords[,1],
  PC2 = pca_coords[,2],
  PC3 = pca_coords[,3],
  cluster = cluster_labels
)
cells_df$text <- paste("Cluster:", as.character(cells_df$cluster))

num_clusters <- length(levels(cluster_labels))
palette <- colorRampPalette(
  c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#EBCC2A","#5BBCD6","#046C9A")
)(num_clusters)
cluster_colors <- palette[as.integer(cells_df$cluster)]

# Simplex vertices
simplex_vertices <- diag(ncol(proj))
vertices_pca <- t(t(simplex_vertices) - pca_res$center) %*% pca_res$rotation[, 1:3]
vertices_df <- data.frame(
  PC1 = vertices_pca[,1],
  PC2 = vertices_pca[,2],
  PC3 = vertices_pca[,3],
  label = paste0("A", 1:k)
)
vertices_df$text <- vertices_df$label

# Cluster slider (All + per cluster)
clusters_with_all <- c("All", levels(cluster_labels))
highlight_archetype_by_cluster <- c("All" = NA)
for (cl in levels(cluster_labels)) {
  idx <- which(cluster_labels == cl)
  alpha_avg <- colMeans(final_model$alphas[idx, , drop = FALSE])
  highlight_archetype_by_cluster[cl] <- which.max(alpha_avg)
}

steps_cluster <- lapply(seq_along(clusters_with_all), function(i) {
  list(
    method = "animate",
    args = list(
      list(clusters_with_all[i]),
      list(mode = "immediate", frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))
    ),
    label = clusters_with_all[i]
  )
})

# Archetype slider
steps_archetype <- lapply(1:k, function(a) {
  cells_k <- which(apply(final_model$alphas, 1, which.max) == a)
  subset_df <- cells_df[cells_k, ]
  subset_colors <- cluster_colors[cells_k]
  
  list(
    method = "restyle",
    args = list(
      list(
        x = list(subset_df$PC1),
        y = list(subset_df$PC2),
        z = list(subset_df$PC3),
        marker = list(size = 3, opacity = 0.6, color = subset_colors)
      ),
      list(0)
    ),
    label = paste0("A", a)
  )
})

# Add "All" option to archetype slider
steps_archetype <- c(
  list(
    list(
      method = "restyle",
      args = list(
        list(
          x = list(cells_df$PC1),
          y = list(cells_df$PC2),
          z = list(cells_df$PC3),
          marker = list(size = 3, opacity = 0.6, color = cluster_colors)
        ),
        list(0)
      ),
      label = "All"
    )
  ),
  steps_archetype
)

# Build figure
fig <- plot_ly()
fig <- add_markers(fig,
                   x = cells_df$PC1,
                   y = cells_df$PC2,
                   z = cells_df$PC3,
                   marker = list(size = 3, opacity = 0.6, color = cluster_colors),
                   text = cells_df$text,
                   hoverinfo = 'text',
                   name = "Cells")

fig <- add_markers(fig,
                   x = vertices_df$PC1,
                   y = vertices_df$PC2,
                   z = vertices_df$PC3,
                   marker = list(size = 8, color = rep("gray", k), symbol = 'diamond'),
                   text = vertices_df$text,
                   hoverinfo = 'text',
                   name = "Archetypes")

# Add simplex edges
edges <- t(combn(1:k, 2))
for (i in seq_len(nrow(edges))) {
  idx1 <- edges[i, 1]; idx2 <- edges[i, 2]
  fig <- add_trace(fig,
                   x = c(vertices_df$PC1[idx1], vertices_df$PC1[idx2]),
                   y = c(vertices_df$PC2[idx1], vertices_df$PC2[idx2]),
                   z = c(vertices_df$PC3[idx1], vertices_df$PC3[idx2]),
                   type = 'scatter3d', mode = 'lines',
                   line = list(color = 'black', width = 2),
                   showlegend = FALSE, hoverinfo = 'none')
}

# Layout with two sliders
fig <- fig %>%
  layout(
    title = 'Archetype Simplex with Cells (colored by cluster)',
    scene = list(
      xaxis = list(title = 'PC1'),
      yaxis = list(title = 'PC2'),
      zaxis = list(title = 'PC3')
    ),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Cluster: "),
        pad = list(t = 60),
        yanchor = "bottom",
        y = -0.1,
        steps = steps_cluster
      ),
      list(
        active = 0,
        currentvalue = list(prefix = "Archetype: "),
        pad = list(t = 20),
        yanchor = "bottom",
        y = -0.3,
        steps = steps_archetype
      )
    )
  ) %>%
  animation_opts(frame = 0, transition = 0, redraw = FALSE)

# Define frames for animations
fig$x$frames <- lapply(clusters_with_all, function(cl) {
  if (cl == "All") {
    df_sub <- cells_df
    colors_sub <- cluster_colors
    archetype_colors <- rep("gray", k)
  } else {
    df_sub <- subset(cells_df, cluster == cl)
    colors_sub <- rep(palette[which(levels(cluster_labels) == cl)], nrow(df_sub))
    hi <- as.integer(highlight_archetype_by_cluster[cl])
    archetype_colors <- rep("gray", k)
    archetype_colors[hi] <- "red"
  }
  
  list(
    name = cl,
    traces = list(0, 1),
    data = list(
      list(
        x = df_sub$PC1, y = df_sub$PC2, z = df_sub$PC3,
        marker = list(size = 3, opacity = 0.6, color = colors_sub),
        text = df_sub$text, hoverinfo = "text",
        mode = "markers", type = "scatter3d", name = "Cells"
      ),
      list(
        x = vertices_df$PC1, y = vertices_df$PC2, z = vertices_df$PC3,
        marker = list(size = 8, color = archetype_colors, symbol = 'diamond'),
        text = vertices_df$text, hoverinfo = "text",
        mode = "markers", type = "scatter3d", name = "Archetypes"
      )
    )
  )
})

fig

### Marker detection analysis by archetype score: save all results to one file
library(Seurat)
library(openxlsx)

# User parameters
k <- 12  # Number of archetypes
default_threshold <- 0.9
min_cells <- 10

# Safety: check that Seurat object exists
if (!exists("seurat_obj")) {
  stop("Seurat object 'seurat_obj' not found in environment.")
}

# Safety: check that required archetype metadata columns exist
arch_cols <- paste0("Archetype_", 1:k)
missing_cols <- setdiff(arch_cols, colnames(seurat_obj@meta.data))
if (length(missing_cols) > 0) {
  stop(paste("Missing archetype metadata columns:", paste(missing_cols, collapse = ", ")))
}

# Initialize Excel workbook
wb <- createWorkbook()

for (i in 1:k) {
  arch_col <- paste0("Archetype_", i)
  
  # Extract archetype scores
  meta_vals <- seurat_obj@meta.data[[arch_col]]
  
  # Dynamically adjust threshold if too few cells
  threshold <- ifelse(sum(meta_vals > default_threshold) < min_cells, 0.8, default_threshold)
  
  cells_near <- rownames(seurat_obj@meta.data)[meta_vals > threshold]
  other_cells <- setdiff(Cells(seurat_obj), cells_near)
  
  # Only run FindMarkers if both groups have enough cells
  if (length(cells_near) >= min_cells && length(other_cells) >= min_cells) {
    # Run marker detection using explicit cell vectors
    markers <- tryCatch({
      FindMarkers(
        object = seurat_obj,
        ident.1 = cells_near,
        ident.2 = other_cells,
        logfc.threshold = 0.25,
        min.pct = 0.1
      )
    }, error = function(e) {
      cat(paste("Error running FindMarkers for", arch_col, ":", e$message, "\n"))
      return(NULL)
    })
    
    if (!is.null(markers) && nrow(markers) > 0) {
      # Add gene names as a column
      markers_out <- cbind(Gene = rownames(markers), markers)
      
      # Add worksheet and write results
      sheet_name <- paste0("Archetype_", i)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, x = markers_out)
      
      cat(paste("Added markers for", sheet_name, "\n"))
    } else {
      cat(paste("No markers found for", arch_col, "\n"))
    }
  } else {
    cat(paste("Skipping", arch_col, "- insufficient cells\n"))
  }
}

# Save results
saveWorkbook(wb, file = "Archetype_Markers.xlsx", overwrite = TRUE)
cat("Saved all archetype markers with gene names included.\n")

# Show unique time points if available
if ("time_point" %in% colnames(seurat_obj@meta.data)) {
  print(unique(seurat_obj$time_point))
} else {
  cat("No 'time_point' column found in metadata.\n")
}
