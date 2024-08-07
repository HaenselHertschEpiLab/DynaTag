# Data Preprocessing
## Demultiplexing via Split Pipe
```bash
# 1 - Demux Sub Libraries and align against hg38_mm39 (although called mm10 in the path)

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-1-cat_S1_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-1-cat_S1_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S1 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-2-cat_S2_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-2-cat_S2_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S2 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-3-cat_S3_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-3-cat_S3_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S3 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-4-cat_S4_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-4-cat_S4_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S4 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-5-cat_S5_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-5-cat_S5_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S5 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-6-cat_S6_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-6-cat_S6_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S6 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-7-cat_S7_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-7-cat_S7_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S7 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-8-cat_S8_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-8-cat_S8_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S8 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
    --mode comb \
    --sublibraries /projects/ag-haensel/Pascal/SCLC/PDX/S1 /projects/ag-haensel/Pascal/SCLC/PDX/S2 /projects/ag-haensel/Pascal/SCLC/PDX/S3 /projects/ag-haensel/Pascal/SCLC/PDX/S4 /projects/ag-haensel/Pascal/SCLC/PDX/S5 /projects/ag-haensel/Pascal/SCLC/PDX/S6 /projects/ag-haensel/Pascal/SCLC/PDX/S7 /projects/ag-haensel/Pascal/SCLC/PDX/S8 \
    --output_dir /projects/ag-haensel/Pascal/SCLC/PDX

# deactivate conda env
conda deactivate
```
# Data Anaylsis via Seurat in R
## Load Packages
```R
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(viridis)
library(reshape2)
library(ComplexHeatmap)
library(DESeq2)
```
## Convenience Functions
```R
rm(list = ls())
# convenience functions
data_path <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration"
fig_path <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration"

SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}
```
## Read in and merge Data
```R
# 1. Reading in Seurat objects
SCLC_PDX_CTRL_1 <- readRDS("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/SCLC_PDX_CTRL_1/DGE_filteredSCLC_PDX_CTRL_1_after_PCA_SCT.RDS")
SCLC_PDX_CTRL_2 <- readRDS("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/SCLC_PDX_CTRL_2/DGE_filteredSCLC_PDX_CTRL_2_after_PCA_SCT.RDS")
SCLC_PDX_CHEM_1 <- readRDS("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/SCLC_PDX_CHEM_1/DGE_filteredSCLC_PDX_CHEM_1_after_PCA_SCT.RDS")
SCLC_PDX_CHEM_2 <- readRDS("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/SCLC_PDX_CHEM_2/DGE_filteredSCLC_PDX_CHEM_2_after_PCA_SCT.RDS")

# Add treatment status
SCLC_PDX_CTRL_1$treatment_status <- "Control"
SCLC_PDX_CTRL_2$treatment_status <- "Control"
SCLC_PDX_CHEM_1$treatment_status <- "Chemotherapy"
SCLC_PDX_CHEM_2$treatment_status <- "Chemotherapy"

# Merge the samples into one Seurat object
merged_object <- merge(SCLC_PDX_CTRL_1, y = list(SCLC_PDX_CTRL_2, SCLC_PDX_CHEM_1, SCLC_PDX_CHEM_2), add.cell.ids = c("SCLC_PDX_CTRL_1", "SCLC_PDX_CTRL_2","SCLC_PDX_CHEM_1", "SCLC_PDX_CHEM_2"))
```
## Dimensionality Reduction
```R
# Re-identify Variable Features
features_SCLC_PDX_CTRL_1 <- VariableFeatures(SCLC_PDX_CTRL_1)
features_SCLC_PDX_CTRL_2 <- VariableFeatures(SCLC_PDX_CTRL_2)
features_SCLC_PDX_CHEM_1 <- VariableFeatures(SCLC_PDX_CHEM_1)
features_SCLC_PDX_CHEM_2 <- VariableFeatures(SCLC_PDX_CHEM_2)

all_features <- Reduce(union, list(features_SCLC_PDX_CTRL_1, features_SCLC_PDX_CTRL_2, features_SCLC_PDX_CHEM_1, features_SCLC_PDX_CHEM_2))
VariableFeatures(merged_object) <- all_features

# Set the Default Assay (if you are working with SCT assay)
DefaultAssay(merged_object) <- "SCT"

# Run PCA on merged data (required before Harmony)
merged_object <- RunPCA(merged_object, features = VariableFeatures(object = merged_object))


# Run ElbowPlot on Harmony-adjusted data
plot <- ElbowPlot(merged_object,ndims = 50)
plot
SaveFigure(plot, "PC_elbow_plot_integrated_SCLC_PDX", width = 8, height = 10)

# Run UMAP on Harmony-adjusted PCA
merged_object <- RunUMAP(merged_object, reduction = "pca", dims = 1:30)

# Find clusters
merged_object <- FindNeighbors(merged_object, reduction = "pca", dims = 1:30)
merged_object <- FindClusters(merged_object, resolution = 0.4)
colours <- DiscretePalette(10, palette = "alphabet", shuffle = FALSE)
DimPlot(merged_object, group.by = "orig.ident", cols = colours) & NoAxes() 
DimPlot(merged_object, group.by = "seurat_clusters", cols = colours, label = T) & NoAxes() & NoLegend()
DimPlot(merged_object, group.by = "treatment_status", cols = colours) & NoAxes()
```
## Assing Cell Types
```R
cell_types <- c("SCLC_0", 
                "SCLC_1", 
                "SCLC_2",
                "SCLC_3",
                "SCLC_4",
                "SCLC_5",
                "Immune Cells",
                "SCLC_7",
                "SCLC_8",
                "Stromal Cells")
names(cell_types) <- levels(merged_object)
merged_object <- RenameIdents(merged_object, cell_types)
CT <- DimPlot(merged_object, cols = colours, label = T) & NoAxes() & NoLegend()
SaveFigure(CT, "DimPlot_CellTypes_SCLC_PDX", width = 8, height = 10)
```
## Remove non-cancerous cells
```R
SCLC_clusters <- c("0","1","2","3","4","5","7","8")
SCLC_object <- subset(x = merged_object, subset = seurat_clusters %in% SCLC_clusters)
DimPlot(SCLC_object, cols = colours, label = T, ) & NoAxes() & NoLegend()
fp <- FeaturePlot(SCLC_object, features = markers, order = F, cols = c("#f5f5f5", "#9b05f2"), keep.scale = "all") & NoAxes()
SaveFigure(fp, "FeaturPlot_TFs_SCLC_PDX", width = 8, height = 10)
```
## Aggregate Counts into Pseudobulk
```R
SCLC_object$orig.ident <- gsub("_", "", SCLC_object$orig.ident)
SCLC_object$info_sample <- paste0(SCLC_object$treatment_status, SCLC_object$orig.ident)

SCLC_CTS <- AggregateExpression(SCLC_object, group.by = "info_sample", assays  = "SCT", return.seurat = F)

AGG_counts <- SCLC_CTS$SCT %>% as.matrix()

colData <- data.frame(samples = colnames(AGG_counts))
colData <- colData %>%
  mutate(condition = ifelse(grepl('Control', samples), 'Control', 'Chemotherapy'))

row.names(colData) <- colData$samples
colData <- colData[, !names(colData) %in% 'samples', drop = FALSE]

colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, ref = "Control")
```
## Set up DESeq2 and run Differential Expression Analysis
```R
dds <- DESeqDataSetFromMatrix(countData = AGG_counts,
                              colData = colData,
                              design = ~condition)

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "condition_Chemotherapy_vs_Control")
res <- na.omit(res)

DESeq2::plotMA(res)
plotDispEsts(dds)
rld <- rlog(dds, blind =T)
```
## PCA
```R
pca_data <- plotPCA(rld, ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
ggsave("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration/PCA.pdf", plot = p, width = 14.72, height = 10.62)
```
## Volcano Plot
```R
res$regulation <- "NS"
res$regulation[res$log2FoldChange > 0.5 & res$padj < 0.1] <- "Up in Chemo" 
res$regulation[res$log2FoldChange < -0.5 & res$padj < 0.1] <- "Down in Chemo" 
up_chemo <- sum(res$regulation == "Up in Chemo")
down_chemo <- sum(res$regulation == "Down in Chemo")

# Adjust the data frame to have negative counts for downregulated genes
regulation_counts <- data.frame(
  Regulation = c("Up in Chemo", "Down in Chemo"),
  Count = c(up_chemo, -down_chemo) # Make the count for "Down in Chemo" negative
)

# Create the bar plot with adjusted counts
ggplot(regulation_counts, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("Up in Chemo" = "blue", "Down in Chemo" = "red")) +
  labs(title = "Number of Genes Upregulated and Downregulated in Chemotherapy",
       x = "Regulation Status",
       y = "Number of Genes (Downregulated shown as negative)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() + # Optional: Flip coordinates to make the bars horizontal
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") # Add a dashed line at y=0

# Create the plot
p <- ggplot(data=as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj), col=regulation, size = abs(log2FoldChange))) + 
  geom_point(alpha=0.75) + 
  theme_minimal() +
  ylab("Significance [-log10(adjusted p-value)]") +
  xlab("Fold Change [log2]") + 
  xlim (-5,5)

# Modify the plot
p2 <- p + 
  scale_color_manual(values= c("black","lightgrey","red")) +
  scale_size_continuous(range = c(0.1, 4)) +  # Adjust the range to make the scale smaller
  theme_classic() +
  labs(color = "Expression", size = "Fold Change Magnitude") +  # Relabeling the legend for color and size
  annotate("text", x = 1, y = 150, label = paste("Up:", sum(res$regulation == "Up in Chemo")), hjust = 0, fontface = "bold", size = 5, colour = "red") +
  annotate("text", x = -1, y = 150, label = paste("Down:", sum(res$regulation  == "Down in Chemo")), hjust = 1, fontface = "bold", size = 5, colour = "black")

p2
ggsave("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration/Volcano.pdf", plot = p2, width = 14.72, height = 10.62)
```
## Heat Map
```R
vsd <- vst(dds)
res_ordered <- res[order(res$padj), ]
top_genes <- rownames(res_ordered)[1:30]
vsd_top <- vsd[top_genes, ]
vsd_top_centered <- t(apply(assay(vsd_top), 1, function(x) x - mean(x)))
colnames(vsd_top_centered) <- c("CHEM_1", "CHEM_2", "CTRL_1", "CTRL_2")
phm <- pheatmap(vsd_top_centered,
         cluster_cols = FALSE, 
         fontsize_row = 8,
         fontsize_col = 8)
# Open a PDF device
pdf("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration/heatmap.pdf")

# Generate the heatmap
phm <- pheatmap(vsd_top_centered,
                cluster_cols = FALSE, 
                fontsize_row = 8,
                fontsize_col = 8)
phm
# Close the PDF device
dev.off()
```
## Pre-ranking for GSEA
```R
res$Rank <- with(res, log2FoldChange * sign(log2FoldChange) * (-log10(padj)))

# Sort genes by rank
res <- res[order(res$log2FoldChange, decreasing = TRUE),]

# Extract ranked gene list and their ranks
ranked_genes <- rownames(res)
ranks <- res$Rank
# Normalize ranks to range between -1 and 1
normalized_ranks <- res$Rank / max(abs(res$Rank))

# Define the number of ranks
num_ranks <- length(ranked_genes)

# Calculate the scaling factor
scaling_factor <- 2 / (num_ranks - 1)

# Calculate the scaled ranks
scaled_ranks <- -1 + (rank(res$Rank) - 1) * scaling_factor

# Write scaled ranks to the .rnk file
write.table(data.frame(Gene = ranked_genes, Rank = scaled_ranks), 
            "SCLC_PDX_Chemo_Control_ranked_genes_scaled.rnk", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
```



























