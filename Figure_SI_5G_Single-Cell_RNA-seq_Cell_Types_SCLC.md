# Preprocessing in RNA-seq_single-cell_SCLC.md
# Plotting of assinged Cell Types in UMAP
```R
merged_object$final_clusters <- merged_object$SCT_snn_res.0.4
merged_object@meta.data[which(merged_object@meta.data$final_clusters == 11),'final_clusters'] <- 3
Idents(merged_object) <- 'final_clusters'
CT <- DimPlot(merged_object, reduction = "umap",raster=FALSE, label=FALSE) 
SaveFigure(CT, "DimPlot_CellTypes_SCLC_PDX", width = 8, height = 10)
```
