# Data Preprocessing
## Download raw BCL data and rename
```bash
bs download run -i 288413178 -o /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/20241203_GP_scDynaTag_mESC_EpiLC

mv 20241203_GP_scDynaTag_mESC_EpiLC BCL_20241203_GP_scDynaTag_mESC_EpiLC
```
## Demultiplex BCL into fastq
```bash
head -30 SampleSheet_v2.csv

[Header],,
FileFormatVersion,2,
RunName,20241203_GP_scDynaTag_mESC_EpiLC,
InstrumentPlatform,NextSeq1k2k,
IndexOrientation,Forward,
,,
[Reads],,
Read1Cycles,51,
Read2Cycles,51,
Index1Cycles,8,
Index2Cycles,8,
,,
[Settings]
SoftwareVersion,4.2.7,
NoLaneSplitting,true,
OverrideCycles,Y51;I8;I8;Y51,
FastqCompressionFormat,gzip,
BarcodeMismatchesIndex1,1,
BarcodeMismatchesIndex2,1,
,,
[Data]
Sample_ID,Index,Index2
R0C0_MYC_EpiLC_1,AACCGGTT,CGTTGGTT,
R0C2_MYC_EpiLC_1,CTGGTCTT,CGTTGGTT,
R0C8_MYC_EpiLC_1,TATACGGA,CGTTGGTT,
R0C21_MYC_EpiLC_1,TACGCGAG,CGTTGGTT,
R0C22_MYC_EpiLC_1,GACTCAAG,CGTTGGTT,
R0C36_MYC_EpiLC_1,CATTCGGT,CGTTGGTT,
R0C41_MYC_EpiLC_1,AAGGTCTG,CGTTGGTT,
R0C62_MYC_EpiLC_1,TAACGCCA,CGTTGGTT,

cat demultiplex_20241203_GP_scDynaTag_mESC_EpiLC.sh

#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --mail-type=END
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Load necessary module
module load bio/bcl2fastq2/2.20.0-GCC-12.2.0

# Define paths for Cell Ranger and inputs
CELLRANGER=/projects/ag-haensel/tools/cellranger-7.0.1/bin/cellranger
RUN_DIR=/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/20241203_GP_scDynaTag_mESC_EpiLC
CSV_FILE=${RUN_DIR}/SampleSheet_v2.csv

# Run Cell Ranger mkfastq
$CELLRANGER mkfastq --id=fastq_20241203_GP_scDynaTag_mESC_EpiLC \
                    --run=$RUN_DIR \
                    --csv=$CSV_FILE
```
## Trim fastq files
```bash
mkdir TRIMMED_DIR

nano trimm_fastq.sh

#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de
# Directory containing the FASTQ files
FASTQ_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path"
# Directory to store trimmed FASTQ files
TRIMMED_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR"
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
## Alignment RHH mm39
```bash
mkdir alignment/sam_mm39

nano alignment_mm39.sh

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de
# Directory containing the trimmed FASTQ files
FASTQ_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR"
# Directory to store sam files
SAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/sam_mm39"
module load bio/Bowtie2/2.5.1-GCC-12.3.0
ref='/projects/ag-haensel/Pascal/genome_files/mm39_bowtie/ref/ref'
# Loop through all R1 FASTQ files
for f1 in $FASTQ_DIR/*_R1_001.trimmed.fastq.gz; do
    # Corresponding R2 file
    f2=$(echo $f1 | sed 's/_R1_001.trimmed.fastq.gz/_R2_001.trimmed.fastq.gz/')
    # Define base_name based on f1
    base_name=$(basename "$f1" _R1_001.trimmed.fastq.gz)
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 -x "$ref" -1 "$f1" -2 "$f2" -S "$SAM_DIR/${base_name}_bowtie2.sam"
done
```
## Alignment RHH mm10
```bash
mkdir alignment/sam_mm10

nano alignment_mm10.sh

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de
# Directory containing the trimmed FASTQ files
FASTQ_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR"
# Directory to store sam files
SAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/sam_mm10"
module load bio/Bowtie2/2.5.1-GCC-12.3.0
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
## Generate BAM files RHH mm39
```bash
mkdir alignment/bam_mm39

nano sam_mm39.to.bam_mm39.sh

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de
# Directory containing the trimmed FASTQ files
SAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/sam_mm39"
# Directory to store sam files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/bam_mm39"
module load bio/SAMtools/1.19.2-GCC-13.2.0
# Loop through all R1 FASTQ files
for f1 in in $SAM_DIR/*_bowtie2.sam; do
    # Define base_name based on f1
    base_name=$(basename "$f1" _bowtie2.sam)
    samtools view -bS -F 0x04 "$f1" > "$BAM_DIR/${base_name}_bowtie2.mapped.bam"
done
```
## Generate BAM files RHH mm10
```bash
mkdir alignment/bam_mm10

nano sam_mm10.to.bam_mm10.sh

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de
# Directory containing the trimmed FASTQ files
SAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/sam_mm10"
# Directory to store sam files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/bam_mm10"
module load bio/SAMtools/1.19.2-GCC-13.2.0
# Loop through all R1 FASTQ files
for f1 in in $SAM_DIR/*_bowtie2.sam; do
    # Define base_name based on f1
    base_name=$(basename "$f1" _bowtie2.sam)
    samtools view -bS -F 0x04 "$f1" > "$BAM_DIR/${base_name}_bowtie2.mapped.bam"
done
```
