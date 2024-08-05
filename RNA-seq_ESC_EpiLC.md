# Data pr-processing of bulk RNA-seq of ESC and EpiLC
## Trimming
```bash
for f in *.fastq.gz; do
sbatch --mem 40GB -J trimming --cpus-per-task 15 --mail-type=FAIL,END --mail-user=phunold@uni-koeln.de --wrap "module load openjdk/11.0.2 && /projects/ag-haensel/tools/bbmap/bbduk.sh in=$f out=${f%%L000_R1_001.fastq.gz}clean.fastq.gz ref=/projects/ag-haensel/tools/bbmap/resources/polyA.fa.gz,/projects/ag-haensel/tools/bbmap/resources/truseq.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20"
done
```
## Alignment via STAR
```bash
#!/bin/bash
    if [ ! -d "STAR_aligned_bam" ]; then
      mkdir STAR_aligned_bam
    fi
    for f in ./*clean.fastq.gz
    do
      module load star/2.7.8a
      STAR --runThreadN 8 \
      --genomeDir /projects/ag-haensel/Michaela/genome_files/mm10_STAR \
      --readFilesIn $f \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.6 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outSAMattributes NH HI NM MD \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix /scratch/phunold/$f \
      --readFilesCommand gunzip -c
    done
```
## Downsampling
```bash
for file in *.bam; do
sbatch --mem 8G -J stat5 --mail-type=FAIL --mail-user=phunold@uni-koeln.de --wrap "module load samtools/1.13 && samtools view -c -F 260 $file >  ${file%%.bam}.stat5"
done

if [ ! -d " downsampled
" ]; then
  mkdir downsampled
fi

for file in *.stat5
do
  stat5_file=`cat $file`
  echo $stat5_file
  if [ "$stat5_file" -gt 8000000 ]
  then
      echo "bigger than 7M" $file
      factor_down=$(awk -v m=$stat5_file 'BEGIN { print 7000000/m }')
      echo $factor_down
      echo "***"
      sbatch --mem 8G --cpus-per-task 8 --wrap "module load samtools/1.13 && samtools view -@ 8 -s $factor_down -b ${file%%.stat5}.bam > ./downsampled/${file%%.stat5}.8M.bam"
  else
    echo " =================> LESS"
    sbatch --mem 8G --wrap "cp ${file%%.stat5}.bam ${file%%.stat5}.less7M.bam"
  fi
  echo "================="
done
```
## Count Reads
```bash
if [ ! -d " counts
" ]; then
  mkdir counts
fi

for bam in *.bam
do
sbatch --mem 16g --mail-type=FAIL,END --mail-user=phunold@uni-koeln.de --wrap "module load use.own && module load pypack/htseq && htseq-count -f bam -m union -s no -t exon -i gene_id $bam /projects/ag-haensel/Michaela/genome_files/mm10_STAR/Mus_musculus.GRCm38.102_edited.gtf > ./counts/${bam%%.merged.bam}.counts.txt"
done
```

# Data Analysis of bulk RNA-seq of ESC and EpiLC

## Load Packages
```R
library(DESeq2)
library(tidyverse)
library(ggplot2)
library("pheatmap")
library("RColorBrewer")
library('org.Mm.eg.db')
library(clusterProfiler)
```
## Read Data and Assemble Count Martix
```R
ESC1 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/ESC-1_S1_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
ESC2 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/ESC-2_S2_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
ESC3 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/ESC-3_S3_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
ESC4 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/ESC-4_S4_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
EpiLC1 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/EpiLC-1_S5_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
EpiLC2 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/EpiLC-2_S6_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
EpiLC3 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/EpiLC-3_S7_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)
EpiLC4 <- read.table('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/EpiLC-4_S8_L001_R1_001.fastq.gzclean.fastq.gzAligned.sortedByCoord.out.7M.bam.counts.txt', header=FALSE, row.names=1, check.names=FALSE)

all(row.names(ESC1) == row.names(ESC2),
    row.names(ESC1) == row.names(ESC3),
    row.names(ESC1) == row.names(ESC4),
    row.names(ESC1) == row.names(EpiLC1),
    row.names(ESC1) == row.names(EpiLC2),
    row.names(ESC1) == row.names(EpiLC3),
    row.names(ESC1) == row.names(EpiLC4)
)

ESC1 <- ESC1
ESC2 <- ESC2[, 1]
ESC3 <- ESC3[, 1]
ESC4 <- ESC4[, 1]
EpiLC1 <- EpiLC1[, 1]
EpiLC2 <- EpiLC2[, 1]
EpiLC3 <- EpiLC3[, 1]
EpiLC4 <- EpiLC4[, 1]

merged_counts_data <- cbind(ESC1, ESC2, ESC3, ESC4, EpiLC1, EpiLC2, EpiLC3, EpiLC4)
merged_counts_data_cleaned <- head(merged_counts_data, n = nrow(merged_counts_data) - 5)
colnames(merged_counts_data_cleaned)[1] <- "ESC1"

colData <- read.delim('/Users/pascalhunold/Desktop/PhD-Documentation/CUT_See/Sequencing/RNAseq/column_data_ESCEpiLC.txt', header=TRUE, row.names=1, check.names=FALSE)

all(colnames(merged_counts_data_cleaned) %in% rownames(colData))
all(colnames(merged_counts_data_cleaned) == rownames(colData))

colData$RA <- factor(colData$RA)

```
## Create DESeq2 Object and Run Differential Expression Analysis
```R
dds <- DESeqDataSetFromMatrix(countData = merged_counts_data_cleaned,
                              colData = colData,
                              design = ~ RA)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$RA <- relevel(dds$RA, ref = "ESC")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resultsNames(dds)

plotMA(res)
```
## Generate Heat Map for top 50 genes
```R
vsd <- vst(dds)
mat <- assay(vsd)[ head(order(res$padj),50), ]
mat <- mat - rowMeans(mat)
#df <- as.data.frame(colData(vsd)[,c("group")])
#pheatmap(mat, annotation_row = df)
pheatmap(mat, cluster_cols = FALSE)

# with gene names
columns(org.Mm.eg.db)
symbols <- rownames(mat)

# use mapIds method to obtain Entrez IDs
ID <- mapIds(org.Mm.eg.db, symbols, 'SYMBOL', 'ENSEMBL')
rownames(mat) <- ID

breaks <- seq(min(mat), max(mat), length.out = 100)
my_palette <- colorRampPalette(c("#000666", "white", "#ffce3f"))(length(breaks)-1)
heatmap <- pheatmap(mat, cluster_cols = FALSE, color = my_palette)
```
## Prepare Data for GSEA
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

















