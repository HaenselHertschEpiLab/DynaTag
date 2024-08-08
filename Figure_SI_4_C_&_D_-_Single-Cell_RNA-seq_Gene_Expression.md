# Preprocessing in RNA-seq_single-cell_SCLC.md
# Plotting of Marker and TF Gene Expression in UMAP
```R
markers <- c("ASCL1", "NEUROD1", "POU2F3", "YAP1", "MYC", "TP53", "TP73", "FOXA1", "NFIB", "NRF1", "SP2")
markers_non_cancer <- c("Ptprc", "Col1a2")

TF <- FeaturePlot(merged_object, features = markers, order = T, cols = c("orange", "darkred")) & NoAxes()
SaveFigure(TF, "DimPlot_TF_SCLC_PDX", width = 8, height = 10)

NC <- FeaturePlot(merged_object, features = markers_non_cancer, order = T, cols = c("orange", "darkred")) & NoAxes()
SaveFigure(NC, "DimPlot_Ptprc_Col1a2_SCLC_PDX", width = 8, height = 10)
```
