# Preprocessing in RNA-seq_single-cell_SCLC.md
#DEA
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
