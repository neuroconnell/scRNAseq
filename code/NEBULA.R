### Here, I deploy NEBULA (NEgative Binomial mixed model Using a Large-sample Approximation, He et al., 2021) to predict differential expression of genes at a cluster resolution, while explicitly modeling treatment and batch as fixed effects.

# single-cell RNA-sequencing datasets are often hierarchically structured, with unequal mean-variance relationships (overdispersion).
# Seurat native functions, when used to identify differentially expressed genes at the cluster-level, treat cells as independent observations leading to artificially small p-values (pseudoreplication).
# As well, replicate level sources of variance are not considered, that is, gene expression may be correlated in some way across all cells derived from the same replicate.

# NEBULA allows us to predict differential expression while explicitly accounting for nested sources of variance, modeled using a negative binomial distribution instead of Poisson distribution. 
# Further, we can predict sources of variance introduced by batch explicitly. Here, "batch" indicates that cDNA from both treatment and control replicates were amplified, libraries were constructed, and replicates were sequenced at the same time.

# We create a negative binomial gamma mixed model where: 
  # 1. replicate-level variation is modeled as a random effect (accounts for overdispersion)
  # 2. treatment and batch are modeled explicitly as fixed effects 
  # 3. random effects are approximated as log-normal (Laplacian large sample approximation)

# We are interested in understanding the impact of a treatment on specific subpopulations of serotonergic neurons, accounting for replicate- and batch-derived variances.


### Defining the NEBULA function
nebula_treat_vs_ctrl_batchasfixed <- function(seurat_obj,
                                              condition_col = "condition",
                                              control_label = "Control",
                                              treatment_label = "Treatment",
                                              cluster_col = "tree.ident") {
  
  Idents(seurat_obj) <- cluster_col
  all_clusters <- sort(as.numeric(as.character(unique(Idents(seurat_obj)))))
  
  nebula_results <- list()
  skipped_clusters <- list()
  
  # Define placeholder sample IDs
  valid_treat_samples <- c("Treatment1", "Treatment2", "Treatment3")
  valid_ctrl_samples  <- c("Control1", "Control2", "Control3")
  all_valid_samples   <- c(valid_treat_samples, valid_ctrl_samples)
  
  # Batch definitions (balanced design)
  batch_map <- c(
    "Treatment1" = "batch_1",
    "Control1"   = "batch_1",
    "Treatment2" = "batch_2",
    "Control2"   = "batch_2",
    "Treatment3" = "batch_3",
    "Control3"   = "batch_3"
  )
  
  for (cl in all_clusters) {
    message("Processing Cluster ", cl)
    
    cluster_cells <- WhichCells(seurat_obj, idents = cl)
    sub_obj <- subset(seurat_obj, cells = cluster_cells)
    meta <- sub_obj@meta.data
    
    meta <- meta %>%
      filter(orig.ident %in% all_valid_samples) %>%
      mutate(
        condition = case_when(
          orig.ident %in% valid_treat_samples ~ treatment_label,
          orig.ident %in% valid_ctrl_samples  ~ control_label,
          TRUE ~ NA_character_
        ),
        condition = factor(condition, levels = c(control_label, treatment_label)),
        replicate = orig.ident,
        batch = factor(batch_map[orig.ident])
      )
    
    # Skip small or unbalanced clusters
    if (nrow(meta) < 50) {
      warning("Skipping cluster ", cl, ": fewer than 50 cells.")
      skipped_clusters[[as.character(cl)]] <- "Fewer than 50 cells"
      next
    }
    if (length(unique(meta$condition)) < 2) {
      warning("Skipping cluster ", cl, ": only one condition group present.")
      skipped_clusters[[as.character(cl)]] <- "Only one condition group"
      next
    }
    if (length(unique(meta$batch)) < 2) {
      warning("Skipping cluster ", cl, ": only one batch group present.")
      skipped_clusters[[as.character(cl)]] <- "Only one batch group"
      next
    }
    
    # Get counts and align metadata
    counts <- GetAssayData(sub_obj, layer = "counts")[, rownames(meta)]
    meta <- meta[rownames(meta), , drop = FALSE]
    
    # Design matrix with condition and batch as fixed effects
    pred <- model.matrix(~ condition + batch, data = meta)
    
    if (any(apply(pred[, -1, drop = FALSE], 2, function(x) length(unique(x))) == 1)) {
      warning("Skipping cluster ", cl, ": zero variance in predictors.")
      skipped_clusters[[as.character(cl)]] <- "Zero-variance predictor"
      next
    }
    
    neb <- tryCatch({
      nebula(
        count     = counts,
        id        = meta$replicate,
        pred      = pred,
        model     = "NBGMM",
        method    = "LN",
        output_re = TRUE
      )
    }, error = function(e) {
      warning("NEBULA failed for cluster ", cl, ": ", e$message)
      skipped_clusters[[as.character(cl)]] <- paste("NEBULA error:", e$message)
      return(NULL)
    })
    
    if (is.null(neb) || is.null(neb$summary)) {
      warning("Skipping cluster ", cl, ": neb$summary missing.")
      skipped_clusters[[as.character(cl)]] <- "NEBULA summary missing"
      next
    }
    
    de <- as.data.frame(neb$summary)
    if (!"gene" %in% colnames(de)) {
      de <- tibble::rownames_to_column(de, "gene")
    }
    
    de <- de %>%
      mutate(FDR = p.adjust(get(paste0("p_condition", treatment_label)), method = "fdr"))
    
    od <- neb$overdispersion %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene")
    
    df <- left_join(de, od, by = "gene") %>%
      mutate(cluster = cl)
    
    nebula_results[[as.character(cl)]] <- df
  }
  
  return(list(results = nebula_results, skipped = skipped_clusters))
}

### Run the model on 
nebula_output <- nebula_treat_vs_ctrl_batchasfixed(
  seurat_obj     = serotonin_obj, 
  condition_col  = "condition",
  control_label  = "Control",
  treatment_label= "Treatment",
  cluster_col    = "tree.ident"
)

### Save results to excel
library(openxlsx)

wb <- createWorkbook()
neb_res <- nebula_output$results
valid_clusters <- names(neb_res)[sapply(neb_res, nrow) > 0]

for (cl in valid_clusters) {
  sheet_name <- paste0("Cluster_", cl)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, neb_res[[cl]])
}

output_file <- "NEBULA_Treat_vs_Control_withBatch.xlsx"
saveWorkbook(wb, file = output_file, overwrite = TRUE)

# Sort Results by batch significance and treatment effect

library(readxl)
library(dplyr)
library(writexl)

input_file  <- "NEBULA_Treat_vs_Control_withBatch.xlsx"
output_file <- "NEBULA_Treat_vs_Control_sorted.xlsx"

sheet_names <- excel_sheets(input_file)
processed_sheets <- list()

for (sheet in sheet_names) {
  df <- read_excel(input_file, sheet = sheet)
  
  treat_col <- grep("^p_conditionTreatment$", colnames(df), value = TRUE)
  batch_cols <- grep("^p_batch", colnames(df), value = TRUE)
  
  if (length(treat_col) != 1 || length(batch_cols) == 0) {
    warning(paste0("Sheet '", sheet, "' missing condition or batch p-values. Skipping."))
    next
  }
  
  df$batch_significant <- apply(df[, batch_cols], 1, function(x) any(x < 0.05, na.rm = TRUE))
  
  df_sorted <- df %>%
    arrange(batch_significant, !!sym(treat_col))
  
  processed_sheets[[sheet]] <- df_sorted
}

write_xlsx(processed_sheets, path = output_file)

#Sort by significant FDR
input_file  <- "NEBULA_Treat_vs_Control_sorted.xlsx"
output_file <- "NEBULA_Treat_vs_Control_FDRsummary.xlsx"

sheet_names <- excel_sheets(input_file)
fdr_summary <- data.frame(sheet = character(),
                          genes_total = integer(),
                          genes_FDRsig = integer(),
                          stringsAsFactors = FALSE)

for (sheet in sheet_names) {
  df <- read_excel(input_file, sheet = sheet)
  fdr_col <- grep("FDR", colnames(df), value = TRUE)
  
  if (length(fdr_col) == 0) {
    warning(paste0("Sheet '", sheet, "' has no FDR column. Skipping."))
    next
  }
  
  fc <- fdr_col[1]
  total_genes <- nrow(df)
  sig_genes   <- sum(df[[fc]] < 0.05, na.rm = TRUE)
  
  fdr_summary <- fdr_summary %>%
    add_row(sheet = sheet, genes_total = total_genes, genes_FDRsig = sig_genes)
  
  message(sprintf("Sheet '%s': %d genes, %d with FDR < 0.05", 
                  sheet, total_genes, sig_genes))
}

write_xlsx(list(FDR_summary = fdr_summary), path = output_file)

                                
