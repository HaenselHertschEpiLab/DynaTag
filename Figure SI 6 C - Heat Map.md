# Preprocessing in RNA-seq_single-cell_SCLC.md
# Heat Map
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
