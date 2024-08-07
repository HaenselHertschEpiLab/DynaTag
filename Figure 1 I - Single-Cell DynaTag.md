# Preprocessing in DynaTag_sc_ESC.md
# Normalisation and Dimentionality Reduction via Seurat in R
## Load Packages
```R
library(ggplot2)
library(Seurat)
```
## Load Data
```R
rm(list = ls())

# Read in the data
counts <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/snDynaTag_ICELL8/Count_matrix/matrix_OCT4.txt")
coordinates <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/snDynaTag_ICELL8/Count_matrix/coordinates_OCT4.txt")

# Set row names for counts
rownames(counts) <- coordinates$V1
```
## Prepare SeuratObject and Normalise
```R
# Create Seurat object
snDynaTag <- CreateSeuratObject(counts = counts)

# Normalize counts
snDynaTag <- NormalizeData(snDynaTag)

# Scale data
snDynaTag <- ScaleData(snDynaTag)
```
## Quality Control
```R
# Calculate the number of peaks per cell
nCount <- colSums(counts)

# Add nCount values as a new column in the Seurat object
snDynaTag[["nCount"]] <- nCount
snDynaTag <- AddMetaData(snDynaTag, metadata = snDynaTag$nCount, col.name = "nCount")

# Plot number of peaks per cell
VlnPlot(snDynaTag, features = c("nCount"), cols = "grey") + 
  NoLegend() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Calculate the number of peaks per cell (non-zero counts)
nPeaks <- apply(counts, 2, function(x) sum(x > 0))

# Add nPeaks values as a new column in the Seurat object
snDynaTag[["nPeaks"]] <- nPeaks
snDynaTag <- AddMetaData(snDynaTag, metadata = snDynaTag$nPeaks, col.name = "nPeaks")

# Plot number of peaks per cell
VlnPlot(snDynaTag, features = c("nPeaks"), cols = "grey") + 
  NoLegend() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Plot nCount and nPeaks on the same plot
VlnPlot(snDynaTag, cols = "grey",
        features = c("nCount", "nPeaks"), ncol = 2)
```
## Filter Data
```R
lower_percentile <- 2
upper_percentile <- 98
lower_cutoff_nPeaks <- quantile(snDynaTag$nPeaks, lower_percentile/100)
upper_cutoff_nPeaks <- quantile(snDynaTag$nPeaks, upper_percentile/100)
snDynaTag <- subset(snDynaTag, subset = nPeaks > lower_cutoff_nPeaks & nPeaks < upper_cutoff_nPeaks)

# Plot nCount and nPeaks after filtering
VlnPlot(snDynaTag, cols = "grey",
        features = c("nCount", "nPeaks"), ncol = 2)
```
## Dimensionality Reduction and Plotting
```R
# Find variable features
snDynaTag <- FindVariableFeatures(snDynaTag)

# Run PCA
snDynaTag <- RunPCA(snDynaTag, features = VariableFeatures(object = snDynaTag))

# Elbow plot
ElbowPlot(snDynaTag, ndims = 40) # Adjust ndims as needed

# Run UMAP
snDynaTag <- RunUMAP(snDynaTag, reduction = "pca", dims = 1:30) # Adjust dims as needed

snDynaTag <- FindNeighbors(snDynaTag, reduction = "pca", dims = 1:30)
snDynaTag <- FindClusters(snDynaTag, resolution = 0.4)
DimPlot(snDynaTag, group.by = "seurat_clusters", label = F, cols = c("blue", "purple", "darkorange")) & NoAxes()# & NoLegend()


# Plot nCount and nPeaks on UMAP
FeaturePlot(snDynaTag, features = "nCount", order = FALSE, cols = c("#f5f5f5", "#9b05f2"), keep.scale = "all") + 
  NoAxes()

FeaturePlot(snDynaTag, features = "nPeaks", order = FALSE, cols = c("#f5f5f5", "#9b05f2"), keep.scale = "all") + 
  NoAxes()
```
