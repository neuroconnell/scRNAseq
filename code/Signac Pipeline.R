### Here, I provide code to analyze and visualize chromatin accessibility using Signac native functions. 
# I load data from the directory, annotate the data, conduct QC and filtering.
# I calculate gene activity and differential accessibility, and then visualize chromatin regions of interest. 

### Signac 
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)   # Change to appropriate genome annotation
library(MACSr)

# Load data from directory
if (!file.exists("matrix.mtx.gz") | !file.exists("barcodes.tsv.gz") | !file.exists("peaks.bed.gz")) {
  stop("One or more required input files (matrix.mtx.gz, barcodes.tsv.gz, peaks.bed.gz) not found.")
}

counts <- Matrix::readMM("matrix.mtx.gz")
barcodes <- readLines("barcodes.tsv.gz")
peaks <- read.table("peaks.bed.gz", sep = "\t", stringsAsFactors = FALSE)
peaknames <- paste(peaks$V1, peaks$V2, peaks$V3, sep = "-")

if (ncol(counts) != length(barcodes)) stop("Mismatch: number of barcodes does not equal matrix columns.")
if (nrow(counts) != nrow(peaks)) stop("Mismatch: number of peaks does not equal matrix rows.")

colnames(counts) <- barcodes
rownames(counts) <- peaknames

fragpath <- "./fragments.tsv.gz"
if (!file.exists(fragpath)) stop("Fragments file not found: ./fragments.tsv.gz")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "GRCm38",
  fragments = fragpath,
  min.cells = 10,
  min.features = 200
)

ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = "seurobj"
)

cat("Seurat object created with", ncol(ATAC), "cells and", nrow(ATAC), "features.\n")

#Supply annotations 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"  # enforce 'chr' style
Annotation(ATAC) <- annotations

# Conduct QC
ATAC <- NucleosomeSignal(ATAC)
ATAC <- TSSEnrichment(ATAC, fast = FALSE)

if (all(c("nucleosome_signal", "TSS.enrichment") %in% colnames(ATAC@meta.data))) {
  DensityScatter(ATAC, x = "nucleosome_signal", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)
  VlnPlot(
    object = ATAC,
    features = c("TSS.enrichment", "nucleosome_signal", "nFeature_peaks"),
    pt.size = 0.1,
    ncol = 3
  )
}

# Filter poor-quality cells
ATAC <- subset(ATAC, subset = TSS.enrichment > 4 & nucleosome_signal > 0.5)
cat("After QC filtering:", ncol(ATAC), "cells remain.\n")

# Normalization and reduction 
ATAC <- RunTFIDF(ATAC)
ATAC <- FindTopFeatures(ATAC, min.cutoff = NULL)
ATAC <- RunSVD(ATAC)
DepthCor(ATAC)

ATAC <- FindNeighbors(ATAC, reduction = "lsi", dims = 2:30)
ATAC <- FindClusters(ATAC, algorithm = 2)
ATAC <- RunUMAP(ATAC, reduction = "lsi", dims = 2:30)

DimPlot(ATAC, label = TRUE) + NoLegend()

# Gene activity analysis
gene.activities <- GeneActivity(ATAC)
ATAC[["RNA"]] <- CreateAssayObject(counts = gene.activities)
ATAC <- NormalizeData(
  ATAC,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(ATAC$nCount_RNA)
)

DefaultAssay(ATAC) <- "RNA"
FeaturePlot(ATAC, features = "Tph2")

# Calculating differential accessibility
DefaultAssay(ATAC) <- "peaks"

if (length(unique(Idents(ATAC))) >= 2) {
  da_peaks <- tryCatch({
    FindMarkers(
      object = ATAC,
      ident.1 = "0",
      ident.2 = "1",
      min.pct = 0.05,
      test.use = "LR"
    )
  }, error = function(e) {
    warning("FindMarkers failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(da_peaks) && nrow(da_peaks) > 0) {
    closest <- ClosestFeature(ATAC, regions = rownames(da_peaks))
    da_peaks$closest_gene <- closest$gene_name
    da_peaks$distance <- closest$distance
    print(head(da_peaks))
  } else {
    cat("No differential peaks found.\n")
  }
} else {
  cat("Not enough clusters to perform DA analysis.\n")
}

# Visualize region of chromatin
roi <- "Gene"
if (roi %in% annotations$gene_name) {
  CoveragePlot(
    object = ATAC,
    region = roi,
    annotation = TRUE,
    peaks = TRUE
  )
} else {
  cat("Gene", roi, "not found in annotations.\n")
}

CoverageBrowser(ATAC, region = 'gene')
