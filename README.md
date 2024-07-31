# DynaTag for efficient profiling of transcription factors in small samples and single cells
This repository contains all code used to analyse and visualise the data from the DynaTag project.

# Code and Data Analysis

## Data Preprocessing (bulk ESC and EpiLC DynaTag)

### Alingment

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

module load bowtie2/2.4.1
ref='/projects/ag-haensel/Pascal/genome_files/mm39_bowtie/ref/ref'

for f1 in *_R1_001.fastq.gz; do
f2=${f1%_R1_001.fastq.gz}_R2_001.fastq.gz
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 -x "$ref" -1 "$f1" -2 "$f2" -S "/scratch/phunold/TFCT/ESC_d2EpiLC_ONSMY/alignment/sam/${f1%%}_bowtie2.sam"
done

for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 8GB --wrap "module load samtools/1.13 && samtools view -bS -F 0x04 "$f1" > /scratch/phunold/TFCT/ESC_d2EpiLC_ONSMY/alignment/bam/${f1%%}_bowtie2.mapped.bam"
done

for f1 in *.bam; do
sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=phunold@uni-koeln.de --wrap "module load samtools/1.13 && samtools view -c -F 260 $f1 >  ${f1%%.bam}.stat5"
done

cat *stat5 > reads.txt

#!/bin/bash

for file in *.bam; do
    if [[ "$file" == *"ATAC"* ]]; then
        desired_reads=15000000
    elif [[  "$file" == *"MYC"* || "$file" == *"YAP1"* ]]; then
        desired_reads=5000000
    elif [[ "$file" == *"OCT4"* ||  ]]; then
        desired_reads=12000000
    elif [[ "$file" == *"SOX2"* ||  ]]; then
        desired_reads=19000000
    elif [[ "$file" == *"H3K27me3"* || "$file" == *"H3K4me3"* || "$file" == *"CTCF"* || "$file" == *"NANOG"* ]]; then
        desired_reads=10000000
    fi
    
    output_file="${file%.*}_${desired_reads}_down_bowtie2.mapped.bam"
    
    sbatch --mem 8G --cpus-per-task 8 --wrap "module load samtools/1.13 && total_reads=\$(samtools view -@ 8 -c -F 260 $file); scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\"); samtools view -@ 8 -bs \"\$scaling_factor\" $file > $output_file"
done

### Peak Calling
### Count Matrix Generation

## Data Preprocessing (single-cell DynaTag)

### Demultiplexing of BCL file
### Alingment
### Peak Calling
### Count and Fragment Matrix Generation

## Data Preprocessing (bulk ESC and EpiLC RNA-seq)

### Alignment 
### Count Matrix Generation

## Data Preprocessing (bulk SCLC PDX DynaTag)

### Aligment
### Peak Calling
### Count Matrix Generation

## Data Preprocessing (single-cell SCLC PDX RNA-seq)

### Demultiplexing and Count Matrix Generation

----

# Figures (Main)

## Figure 1
### Figure 1 B - Correlation Plot

### Figure 1 C - FRiP Score

### Figure 1 D - IGV Snapshot

### Figure 1 E - MEME-ChIP Motif Prediction

### Figure 1 F - FRiP Score

### Figure 1 G - Differential Binding Analysis (edgeR)

### Figure 1 H - TOBIAS Motif Prediction

### Figure 1 I - Single-Cell DynaTag

## Figure 2
### Figure 2 B - Differential Binding Analysis (DiffBind)

### Figure 2 C - PlotProfile

### Figure 2 D - Bubble Plot


# Figures (Supplementary Information)

## Figure SI 1
### Figure SI 1 B - Correlation Plot

### Figure SI 1 C - FRiP Score

### Figure SI 1 D - Target Gene Enrichment

## Figure SI 2 
### Figure SI 2 A - Differential Binding Analysis (edgeR)

### Figure SI 2 C - GSEA

## Figure SI 3 
### Figure SI 3 A - FRiP Score

### Figure SI 3 B - Target Gene Enrichment

### Figure SI 3 C - Single-Cell DynaTag

## Figure SI 4
### Figure SI 4 - Differential Binding Analysis (DiffBind)

## Figure SI 5
### Figure SI 5 A - Single-Cell RNA-seq QC

### Figure SI 5 B - Single-Cell RNA-seq Marker Gene Expression

### Figure SI 5 C - Single-Cell RNA-seq TF Expression

### Figure SI 5 D - Single-Cell RNA-seq Cell Types

## Figure SI 6 
### Figure SI 6 A - PCA

### Figure SI 6 B - Differential Expression Analysis (DESeq2)

### Figure SI 6 C - Heat Map

### Figure SI 6 D - GSEA


















