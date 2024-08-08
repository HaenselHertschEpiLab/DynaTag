# Preprocessing in RNA-seq_single-cell_SCLC.md
# Plotting of assinged Cell Types in UMAP
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
