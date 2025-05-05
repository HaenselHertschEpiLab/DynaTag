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
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(viridis)
library(reshape2)
library(ComplexHeatmap)
library(DESeq2)
library(stringr)
library(purrr)
library(scclusteval)
library(tidyr)
library(parallel)
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
## Pre-processing
```R
# Predict cell cycle 
merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])
merged_object <- NormalizeData(merged_object)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_object<-CellCycleScoring(merged_object, g2m.features=g2m.genes[g2m.genes %in% rownames(merged_object)], s.features=s.genes[s.genes %in% rownames(merged_object)], set.ident = FALSE)
merged_object@meta.data$CC.Difference <- merged_object@meta.data$S.Score - merged_object@meta.data$G2M.Score

# Split the Seurat object by sample for preprocessing
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f = merged_object$orig.ident)

# QC plots
pdf('VlnPlot.data_nGene_filtered.pdf')
VlnPlot(merged_object, features="nFeature_RNA", pt.size = 0.3, group.by= "orig.ident")
dev.off()
pdf('VlnPlot.data_nUMI_filtered.pdf')
VlnPlot(merged_object, c("nCount_RNA"), pt.size = 0.3, group.by= "orig.ident")
dev.off()
merged_object <- PercentageFeatureSet(merged_object, pattern = "^MT-", col.name = "percent.mito")
pdf('ViolinPlots.expr_percMito_filtered.pdf')
VlnPlot(merged_object, c("percent.mito"),  pt.size = 0.3, group.by= "orig.ident")
dev.off()
merged_object <- PercentageFeatureSet(merged_object, pattern = "^RPL", col.name = "percent.ribo")
pdf('ViolinPlots.expr_percRibo_filtered.pdf')
VlnPlot(merged_object, c("percent.ribo"),  pt.size = 0.3, group.by= "orig.ident")
dev.off()

# Filter cells 
merged_object <- subset(merged_object, subset= nFeature_RNA > 500 & nCount_RNA > 1000 & nCount_RNA < 100000 )
merged_object <- subset(merged_object, subset= percent.mito < 10)

```
## Integration and clustering
```R
# Normalization and PCA
merged_object <- SCTransform(merged_object, vars.to.regress = c("percent.mito"))
merged_object <- RunPCA(merged_object, pcs.compute=50, assay = "SCT")
# Evaluate relevant PC
pdf('ElbowPlot.pdf')
ElbowPlot(merged_object, ndims = 50)
dev.off()

# RPCA integration ans clustering 
merged_object <- IntegrateLayers(object = merged_object,method = RPCAIntegration,normalization.method = "SCT",new.reduction = "integrated.rpca",verbose = F)
merged_object <- FindNeighbors(merged_object, dims = 1:30, reduction = "integrated.rpca")
merged_object <- FindClusters(merged_object, resolution = c(0.2,0.3,0.4,0.5))
merged_object <- RunUMAP(merged_object, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap")

pdf('umap_plot_clusters_rpca.pdf')
for ( i in colnames(merged_object@meta.data)[grep('SCT_snn',colnames(merged_object@meta.data))]){
  Idents(merged_object) <- i
  print(DimPlot(merged_object, reduction = "umap",raster=FALSE, label=TRUE) + ggtitle(paste(i)))
}
dev.off()

```
## Cluster stability analysis
```R
# Cluster stability analysis (for each resolution)
obj<-merged_object
obj<-AddMetaData(obj, merged_object@meta.data[,grep("SCT_snn_res", colnames(merged_object@meta.data))])
resolutions <- colnames(merged_object@meta.data)[grep("SCT_snn_res", colnames(merged_object@meta.data))]
# This function re-compute integration and clustering after each subsampling
performSubsamplingAndClustering <- function(obj, resolutions, nPCs, perc.sub=0.8, n_subsampling=20){
  subsample_idents<-Reduce( "bind_rows", lapply(resolutions, function(x){
    res<-str_remove(x,"SCT_snn_res.")
    out<-mclapply(c(1:n_subsampling), function(y){
      rand_test <- RandomSubsetData(obj, rate=perc.sub)
      rand_test_rpca <- IntegrateLayers(object = rand_test,method = RPCAIntegration,normalization.method = "SCT",new.reduction = "integrated.rpca",verbose = F)
      rand_test_rpca <- FindNeighbors(rand_test_rpca, dims = 1:nPCs, verbose = FALSE,reduction = "integrated.rpca")
      eval(parse(text=paste0("rand_test_rpca <- FindClusters(rand_test_rpca, resolution = ",res,", verbose = FALSE)")))
      recluster_ident<-rand_test_rpca@meta.data[,x]
      names(recluster_ident) <- rownames(rand_test_rpca@meta.data)
      original_ident<-obj@meta.data[rownames(rand_test_rpca@meta.data),x]
      names(original_ident) <- rownames(rand_test_rpca@meta.data)
      return(list(recluster_ident,original_ident))
    },mc.cores = 6)
    recluster_ident <-list()
    original_ident <-list()
    for(i in 1:n_subsampling){
      recluster_ident[[i]]=out[[i]][[1]]
      original_ident[[i]]=out[[i]][[2]]
    }
    data<-tibble(resolution=res, recluster_ident=recluster_ident, original_ident=original_ident,round=as.character(c(0:(n_subsampling-1))))
    return(data)
  }))
  return(subsample_idents)
}
subsample_idents <- performSubsamplingAndClustering(obj, resolutions,nPCs=30)

# Find most stable clustering
subsample_idents_list<- subsample_idents %>% group_by(resolution) %>%  nest()
stable_cl_res0.2<-AssignStableCluster(subsample_idents_list$data[[1]]$original_ident,subsample_idents_list$data[[1]]$recluster_ident,method = "jaccard_median",jaccard_cutoff = 0.75)
stable_cl_res0.3<-AssignStableCluster(subsample_idents_list$data[[2]]$original_ident,subsample_idents_list$data[[2]]$recluster_ident,method = "jaccard_median",jaccard_cutoff = 0.75)
stable_cl_res0.4<-AssignStableCluster(subsample_idents_list$data[[3]]$original_ident,subsample_idents_list$data[[3]]$recluster_ident,method = "jaccard_median",jaccard_cutoff = 0.75)
stable_cl_res0.5<-AssignStableCluster(subsample_idents_list$data[[4]]$original_ident,subsample_idents_list$data[[4]]$recluster_ident,method = "jaccard_median",jaccard_cutoff = 0.75)

barplot_data<-data.frame(fraction = c(sum(stable_cl_res0.2$stable_cluster)/length(stable_cl_res0.2$stable_cluster),sum(stable_cl_res0.3$stable_cluster)/length(stable_cl_res0.3$stable_cluster),sum(stable_cl_res0.4$stable_cluster)/length(stable_cl_res0.4$stable_cluster),sum(stable_cl_res0.5$stable_cluster)/length(stable_cl_res0.5$stable_cluster)), resolution=c("res0.2","res0.3","res0.4","res0.5"))
pdf('Barplot_fraction_of_stable_clusters.pdf')
ggplot(data=barplot_data, aes(x=resolution, y=fraction)) + geom_bar(stat="identity", fill="steelblue") + ylim(0,1) + ylab("Fraction of stable clusters") + xlab("Resolution") + theme_minimal() + geom_text(aes(label=c("n = 8","n = 9","n = 11","n = 11")), vjust=1.6, color="white", size=5)
dev.off()

JI_res0.2<-AssignHighestJaccard(subsample_idents_list$data[[1]]$data[[1]]$original_ident, subsample_idents_list$data[[1]]$data[[1]]$recluster_ident)
JI_res0.3<-AssignHighestJaccard(subsample_idents_list$data[[2]]$data[[1]]$original_ident, subsample_idents_list$data[[2]]$data[[1]]$recluster_ident)
JI_res0.4<-AssignHighestJaccard(subsample_idents_list$data[[3]]$data[[1]]$original_ident, subsample_idents_list$data[[3]]$data[[1]]$recluster_ident)
JI_res0.5<-AssignHighestJaccard(subsample_idents_list$data[[4]]$data[[1]]$original_ident, subsample_idents_list$data[[4]]$data[[1]]$recluster_ident)
data_boxplot <- data.frame(median_JI = c(apply(JI_res0.2,2,median),apply(JI_res0.3,2,median),apply(JI_res0.4,2,median),apply(JI_res0.5,2,median)), Resolution=c(rep("res0.2",length(apply(JI_res0.2,2,median))),rep("res0.3",length(apply(JI_res0.3,2,median))),rep("res0.4",length(apply(JI_res0.4,2,median))),rep("res0.5",length(apply(JI_res0.5,2,median)))))
pdf('Boxplot_median_JI.pdf')
ggplot(data_boxplot, aes(x=Resolution, y=median_JI)) + geom_boxplot(outliers = FALSE, color="steelblue") + geom_jitter(shape=16, position=position_jitter(0.2)) + ylab("Clusters median JI") + geom_hline(yintercept = 0.75,colour="red3",size=0.8,linetype=8 ) + theme_minimal()
dev.off()

pdf('Jaccard_plot_res0.4_stable_cluster.pdf')
JaccardRainCloudPlot(idents1 = subsample_idents_list$data[[3]]$original_ident,
                     idents2 = subsample_idents_list$data[[3]]$recluster_ident) + 
  geom_hline(yintercept = c(0.75), linetype = 2) + geom_hline(yintercept = c(0.75), linetype = 2,color="red2") +
  xlab("clusters res=0.4") 
dev.off()

```
## Cluster identity
```R
# After evaluating cluster stability we keep res 0.4 (merging clusters 3 and 11, because splitting is not stable)
merged_object$final_clusters <- merged_object$SCT_snn_res.0.4
merged_object@meta.data[which(merged_object@meta.data$final_clusters == 11),'final_clusters'] <- 3
pdf('umap_plot_final_clusters_rpca.pdf')
Idents(merged_object) <- 'final_clusters'
DimPlot(merged_object, reduction = "umap",raster=FALSE, label=FALSE) 
dev.off()

## Compute cluster markers
merged_object<- PrepSCTFindMarkers(merged_object, assay = "SCT", verbose = TRUE)
Idents(merged_object) <- 'final_clusters'
DefaultAssay(merged_object) <- 'SCT'
MarkersAll <- Reduce("rbind",lapply(levels(Idents(merged_object)), function(cluster) {
  Markers <- FindMarkers(merged_object, ident.1 = cluster, ident.2 = NULL, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 1, pseudocount.use = 0.2,assay = "SCT")
  Markers <- Markers[which(Markers$p_val_adj < 0.01),]
  Markers <- Markers[order(Markers$avg_log2FC, decreasing = TRUE),]
  Markers$gene <- rownames(Markers)
  Markers$cluster <- rep(paste0(cluster),nrow(Markers))
  print(nrow(Markers))
  return(Markers)
}))

# Prepare dot plots for cluster markers
mito_genes <- grep(pattern = "^MT", x = rownames(merged_object), value = TRUE)
markers<-MarkersAll[!(MarkersAll$gene %in% mito_genes),]
markers$logFC <-markers$avg_log2FC
markers<-markers %>% group_by(cluster) %>% top_n(n = 10)
markers<- markers[order(as.numeric(markers$cluster), decreasing = FALSE),]
p <- DotPlot(object = merged_object, assay="SCT", features=unique(markers$gene), scale=TRUE,idents = c(0:10),dot.scale = 9) + ylab("Clusters") + scale_colour_gradient2(low = "blue2", mid = "gray90", high = "red2", midpoint=0) + scale_size_continuous(limits = c(0, 100))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
pdf('DotPlot_final_clusters_top10_markers_rank_avg_logFC.pdf', width=20)
print(p)
dev.off()

# UMAP plots with gene expresssion
pdf('umap_plot_Ptprc_expr.pdf')
print(FeaturePlot(merged_object, features = "Ptprc", pt.size = 0.4,slot = "data", order=TRUE) + scale_colour_gradientn(colours = c('grey90',"orange","red4")))
dev.off()

pdf('umap_plot_Col1a2_expr.pdf')
print(FeaturePlot(merged_object, features = "Col1a2", pt.size = 0.4,slot = "data", order=TRUE) + scale_colour_gradientn(colours = c('grey90',"orange","red4")))
dev.off()

## Plot for TFs
TF_names<-c("ASCL1","YAP1","TP53","TP73","FOXA1","NFIB","NRF1","SP2")
pdf('umap_plot_TFs_expr.pdf')
for (gene in TF_names){
  print(FeaturePlot(merged_object, features = gene, pt.size = 0.4,slot = "data", order=TRUE) + scale_colour_gradientn(colours = c('grey90',"orange","red4")))
}
dev.off()
```
## Remove non-cancerous cells
```R
SCLC_clusters <- c("0","1","2","3","4","5","7","8")
SCLC_object <- subset(x = merged_object, subset = final_clusters %in% SCLC_clusters)
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



























