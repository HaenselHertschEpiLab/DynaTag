# Preprocessing in RNA-seq_ESC_EpiLC.md
# Prepare Data for GSEA
```R
# keep only necessary columns
padj_log <- as.data.frame(cbind(res$padj, res$log2FoldChange))
#rename columns
colnames(padj_log) <- c("padj","log2FoldChange")
# keep gene names as rownames
rownames(padj_log) <- rownames(res)
# create a colum with the information if up or downregualted
padj_log$sign <- ifelse(padj_log$log2FoldChange >0, "1", "-1")
# replace na with 0
padj_log[is.na(padj_log)] <- 0
padj_log$sign <- as.numeric(padj_log$sign)

# create rank column:
# pvalue times sign (log2FoldChange up or down from before)
padj_log$rank <- (padj_log$padj*padj_log$sign)

#order by rank
padj_log_ord <- padj_log[order(-padj_log$rank),]

# make sure it's ranked correctly
pos <- padj_log_ord[padj_log_ord$rank>0,]
null <- padj_log_ord[padj_log_ord$rank==0,]
neg <- padj_log_ord[padj_log_ord$rank<0,]

pos$rank2 <- 1-pos$rank
neg$rank2 <- -(1+neg$rank)
null$rank2 <- null$rank

# combine to final ranked list
rank_final <- rbind(pos,null,neg)
rank_final$ID <- rownames(rank_final)
rank_final <- rank_final[order(-rank_final$rank2),]

ranks <- rank_final[,c(6,5)]

#save
write.table(ranks, "/Users/pascalhunold/Desktop/PhD-Documentation/CUTSee/Sequencing/RNAseq/ESCvsEpiLC_GSEA_ranked_genes.rnk", quote = F, row.names = F, col.names = F, sep="\t")
```
