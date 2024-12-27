# Data Preprocessing
## Download raw BCL data
```bash
bs download run -i 272592385 -o /scratch/phunold/20231231_icell8cx_ESC
```
## Demultiplex BCL into fastq
```bash
#!/bin/bash -l
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --mail-type=END
#SBATCH --mail-user=phunold@uni-koeln.de
module load use.own
module load module load bcl2fastq/2.20.0_rh7
/projects/ag-haensel/tools/cellranger-7.0.1/bin/cellranger mkfastq --id=20231231_icell8cx_ESC \
                     --run=/scratch/phunold/20231231_icell8cx_ESC \
                     --csv=/scratch/phunold/20231231_icell8cx_ESC/SampleSheet_UDP.csv
```
## RHH Demultiplex BCL into fastq
```bash
head SampleSheet_UDP_RHH.csv

[Header],,
FileFormatVersion,2,
RunName,PH_20231231_snDynaTag_ICELL8CX_ESC_UDP,
InstrumentPlatform,NextSeq1k2k,
IndexOrientation,Forward,
,,
[Reads],,
Read1Cycles,50,
Read2Cycles,50,
Index1Cycles,10,
Index2Cycles,10,
,,
[Settings]
SoftwareVersion,3.10.4,
NoLaneSplitting,true,
OverrideCycles,Y50;I10;I10;Y50,
FastqCompressionFormat,gzip,
BarcodeMismatchesIndex1,1,
BarcodeMismatchesIndex2,1,
,,
[Data]
Sample_ID,Index,Index2
R0C2_OCT4_ESC,GACTGAGTAG,CGCTCCACGA,
R0C22_OCT4_ESC,ACCGGCCGTA,CGCTCCACGA,
R0C24_OCT4_ESC,CTGCGAGCCA,CGCTCCACGA,
#... --> Lines are missing, for each cell that was evaluated as single nucleus by the ICELL8 software, a unique line #in the sample sheet was created indicating the Well poition (RXCX) and transcripton factor targeted. 

cat demultiplex.sh

#!/bin/bash
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --mail-type=END
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Load necessary module
module load bio/bcl2fastq2/2.20.0-GCC-12.2.0

# Define paths for Cell Ranger and inputs
CELLRANGER=/projects/ag-haensel/tools/cellranger-7.0.1/bin/cellranger
RUN_DIR=/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/Pascal_snDynaTag/snDynaTag_ICELL8cx_BCL
CSV_FILE=${RUN_DIR}/SampleSheet_UDP_RHH.csv

# Run Cell Ranger mkfastq
$CELLRANGER mkfastq --id=20231231_icell8cx_ESC \
                    --run=$RUN_DIR \
                    --csv=$CSV_FILE
```
## Trim fastq files
```bash
mkdir TRIMMED_DIR

#!/bin/bash
    #SBATCH --time=6:00:00
    #SBATCH --mem=32gb
    #SBATCH --mail-type=ALL
    #SBATCH --mail-user=phunold@uni-koeln.de
# Directory containing the FASTQ files
FASTQ_DIR="/scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path"
# Directory to store trimmed FASTQ files
TRIMMED_DIR="/scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR"
module load use.own && module load pypack/cutadapt
# Loop through all R1 FASTQ files
for f1 in $FASTQ_DIR/*_R1_001.fastq.gz; do
    # Corresponding R2 file
    f2=$(echo $f1 | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
    # Extract the basename for naming
    base_name=$(basename $f1 _R1_001.fastq.gz)
    # Run cutadapt directly
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 \
             -o $TRIMMED_DIR/${base_name}_R1_001.trimmed.fastq.gz \
             -p $TRIMMED_DIR/${base_name}_R2_001.trimmed.fastq.gz \
             $f1 $f2
done
```
## Trim fastq files RHH 27122024
```bash
mkdir TRIMMED_DIR

nano trimm_fastq.sh

#!/bin/bash
    #SBATCH --time=6:00:00
    #SBATCH --mem=8gb
# Directory containing the FASTQ files
FASTQ_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/20231231_icell8cx_ESC/outs/fastq_path"
# Directory to store trimmed FASTQ files
TRIMMED_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR"
module load bio/cutadapt/4.4-GCCcore-12.2.0
# Loop through all R1 FASTQ files
for f1 in $FASTQ_DIR/*_R1_001.fastq.gz; do
    # Corresponding R2 file
    f2=$(echo $f1 | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
    # Extract the basename for naming
    base_name=$(basename $f1 _R1_001.fastq.gz)
    # Run cutadapt directly
    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 \
             -o $TRIMMED_DIR/${base_name}_R1_001.trimmed.fastq.gz \
             -p $TRIMMED_DIR/${base_name}_R2_001.trimmed.fastq.gz \
             $f1 $f2
done
```
## Alignment
```bash
mkdir /alignment/sam

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
# Directory containing the trimmed FASTQ files
FASTQ_DIR="/scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR"
# Directory to store sam files
SAM_DIR="/scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR/alignment/sam"
module load bowtie2/2.4.1
ref='/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/ref/ref'
# Loop through all R1 FASTQ files
for f1 in $FASTQ_DIR/*_R1_001.trimmed.fastq.gz; do
    # Corresponding R2 file
    f2=$(echo $f1 | sed 's/_R1_001.trimmed.fastq.gz/_R2_001.trimmed.fastq.gz/')
    # Define base_name based on f1
    base_name=$(basename "$f1" _R1_001.trimmed.fastq.gz)
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 -x "$ref" -1 "$f1" -2 "$f2" -S "$SAM_DIR/${base_name}_bowtie2.sam"
done
```
## Generate BAM files
```bash
mkdir /alignment/bam

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
# Directory containing the trimmed FASTQ files
SAM_DIR="/scratch/phunold/icell8cx_20231212/mkfastq/snDynaTag_ICELL8cx_ESC_EpiLC_20231212/TRIMMED_DIR/alignment/sam"
# Directory to store sam files
BAM_DIR="/scratch/phunold/icell8cx_20231212/mkfastq/snDynaTag_ICELL8cx_ESC_EpiLC_20231212/TRIMMED_DIR/alignment/bam"
module load samtools/1.13
# Loop through all R1 FASTQ files
for f1 in in $SAM_DIR/*_bowtie2.sam; do
    # Define base_name based on f1
    base_name=$(basename "$f1" _bowtie2.sam)
    samtools view -bS -F 0x04 "$f1" > "$BAM_DIR/${base_name}_bowtie2.mapped.bam"
done
```
## Merge single-cell BAM files into aggregated BAM file
```bash
## adjust row and col parameters according to sample coordinates on ICELL8 chip

#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

# Directory containing BAM files
BAM_DIR="/scratch/phunold/icell8cx_20231212/mkfastq/snDynaTag_ICELL8cx_ESC_EpiLC_20231212/TRIMMED_DIR/alignment/bam"

# Output directory for merged BAM file
OUTPUT_DIR="/scratch/phunold/icell8cx_20231212/mkfastq/snDynaTag_ICELL8cx_ESC_EpiLC_20231212/TRIMMED_DIR/alignment/aggregated_BAM"
mkdir -p $OUTPUT_DIR

# Output file name
OUTPUT_FILE="${OUTPUT_DIR}/merged_sample_EpiLC_NANOG.bam"

# Create an empty array to store file names
declare -a bam_files_to_merge

# Loop through the specified range of rows and columns
for row in {18..36}; do
    for col in {36..71}; do
        # Construct the file pattern
        file_pattern="R${row}C${col}*_bowtie2.mapped.bam"

        # Add matching files to the array
        for file in ${BAM_DIR}/${file_pattern}; do
            if [ -f "$file" ]; then
                bam_files_to_merge+=("$file")
            fi
        done
    done
done

# Merge the files using samtools
module load samtools/1.13
samtools merge -@ 8 $OUTPUT_FILE "${bam_files_to_merge[@]}"
```
## Peak Calling via MACS2 from aggregated BAM files
```bash
for f in *sorted.bam; do
sbatch -J MACS2 --mem 32GB --wrap "module load use.own && module load pypack/macs2 && macs2 callpeak -t $f -f BAMPE -g mm --keep-dup all -n /scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_BAM/peaks/$f --nomodel --extsize 55 -B --SPMR"
done

module load bedtools/2.29.2
for f in *_peaks.narrowPeak; do
  awk '{print $1"\t"$2"\t"$3}' $f | sortBed -i - | mergeBed -i - > ${f%%.bed}_real.bed
done
```
## Generate master peak file
```bash
module load bedtools/2.29.2
touch master_peak.bed
for f in *_real.bed; do
sortBed -i $f | mergeBed -i - >> merged_peak.bed
done
sortBed -i merged_peak.bed | mergeBed -i - > master_peak.bed
awk '($3 - $2) > 59' master_peak.bed > master_peak_60nt.bed
```
## Generate Count Matrix per Cell
```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module load bedtools/2.29.2
for f1 in R*.bam; do
out_name="$(basename "${f1}" .bam).bedgraph"
bedtools coverage -a /scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_BAM/peaks/master_peak_60nt.bed -b "${f1}" -counts > /scratch/phunold/20231231_icell8cx_ESC/20231231_icell8cx_ESC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_BAM/peaks/bedgraph/"${out_name}"
done

#!/bin/bash
for f1 in *.bedgraph; do
    base_name=$(basename "$f1" .bedgraph)
    txt_file="${base_name}.txt"
    awk '{print "chr"$1":"$2"-"$3 "\t" $4}' "$f1" > "$txt_file"
done

#!/bin/bash

txt_files=(R*.txt)
prev_first_column=""
first_file=true

for txt_file in "${txt_files[@]}"; do
current_first_column=$(cut -f1 "$txt_file")
if [ "$first_file" = true ]; then
first_file=false
prev_first_column="$current_first_column"
continue
fi
if [ "$current_first_column" != "$prev_first_column" ]; then
echo "Mismatch found in $txt_file. The first column is different from the previous file."
fi 
prev_first_column="$current_first_column"
done
echo "Completed checking for mismatches."
```
# Normalisation and Dimentionality Reduction via Seurat in R
## Load Packages
```R
library(ggplot2)
library(Seurat)
```
## Load Data
```R
rm(list = ls())

# Read in the data
counts <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/snDynaTag_ICELL8/Count_matrix/matrix_OCT4.txt")
coordinates <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/snDynaTag_ICELL8/Count_matrix/coordinates_OCT4.txt")

# Set row names for counts
rownames(counts) <- coordinates$V1
```
## Prepare SeuratObject and Normalise
```R
# Create Seurat object
snDynaTag <- CreateSeuratObject(counts = counts)

# Normalize counts
snDynaTag <- NormalizeData(snDynaTag)

# Scale data
snDynaTag <- ScaleData(snDynaTag)
```
## Quality Control
```R
# Calculate the number of peaks per cell
nCount <- colSums(counts)

# Add nCount values as a new column in the Seurat object
snDynaTag[["nCount"]] <- nCount
snDynaTag <- AddMetaData(snDynaTag, metadata = snDynaTag$nCount, col.name = "nCount")

# Plot number of peaks per cell
VlnPlot(snDynaTag, features = c("nCount"), cols = "grey") + 
  NoLegend() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Calculate the number of peaks per cell (non-zero counts)
nPeaks <- apply(counts, 2, function(x) sum(x > 0))

# Add nPeaks values as a new column in the Seurat object
snDynaTag[["nPeaks"]] <- nPeaks
snDynaTag <- AddMetaData(snDynaTag, metadata = snDynaTag$nPeaks, col.name = "nPeaks")

# Plot number of peaks per cell
VlnPlot(snDynaTag, features = c("nPeaks"), cols = "grey") + 
  NoLegend() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Plot nCount and nPeaks on the same plot
VlnPlot(snDynaTag, cols = "grey",
        features = c("nCount", "nPeaks"), ncol = 2)
```
## Filter Data
```R
lower_percentile <- 2
upper_percentile <- 98
lower_cutoff_nPeaks <- quantile(snDynaTag$nPeaks, lower_percentile/100)
upper_cutoff_nPeaks <- quantile(snDynaTag$nPeaks, upper_percentile/100)
snDynaTag <- subset(snDynaTag, subset = nPeaks > lower_cutoff_nPeaks & nPeaks < upper_cutoff_nPeaks)

# Plot nCount and nPeaks after filtering
VlnPlot(snDynaTag, cols = "grey",
        features = c("nCount", "nPeaks"), ncol = 2)
```
## Dimensionality Reduction and Plotting
```R
# Find variable features
snDynaTag <- FindVariableFeatures(snDynaTag)

# Run PCA
snDynaTag <- RunPCA(snDynaTag, features = VariableFeatures(object = snDynaTag))

# Elbow plot
ElbowPlot(snDynaTag, ndims = 40) # Adjust ndims as needed

# Run UMAP
snDynaTag <- RunUMAP(snDynaTag, reduction = "pca", dims = 1:30) # Adjust dims as needed

snDynaTag <- FindNeighbors(snDynaTag, reduction = "pca", dims = 1:30)
snDynaTag <- FindClusters(snDynaTag, resolution = 0.4)
DimPlot(snDynaTag, group.by = "seurat_clusters", label = F, cols = c("blue", "purple", "darkorange")) & NoAxes()# & NoLegend()


# Plot nCount and nPeaks on UMAP
FeaturePlot(snDynaTag, features = "nCount", order = FALSE, cols = c("#f5f5f5", "#9b05f2"), keep.scale = "all") + 
  NoAxes()

FeaturePlot(snDynaTag, features = "nPeaks", order = FALSE, cols = c("#f5f5f5", "#9b05f2"), keep.scale = "all") + 
  NoAxes()
```


























































