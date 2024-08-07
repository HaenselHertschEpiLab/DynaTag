# Preprocessing in RNA-seq_single-cell_SCLC.md
# PCA
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
