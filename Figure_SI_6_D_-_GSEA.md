# Preprocessing in RNA-seq_single-cell_SCLC.md
# Pre-ranking for GSEA
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
