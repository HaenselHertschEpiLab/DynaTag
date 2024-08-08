# Preprocessing in RNA-seq_single-cell_SCLC.md
# Plot QC Metrics as Violins in R
```R
features <- c("nFeature_RNA", "nCount_RNA")
plot <- VlnPlot(merged_object, features = features, pt.size = 0.10, ncol = 2, cols = c("grey","grey","grey","grey"))
SaveFigure(plot, "QC_nFeature_nCount", width = 8, height = 10)

```
