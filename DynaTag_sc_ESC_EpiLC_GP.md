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
## Merge single-cell BAM files into aggregated BAM file mm39
```bash
nano aggregated_bam_mm39.sh

#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Directory containing BAM files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/bam_mm39"

# Output directory for merged BAM files
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm39"
mkdir -p "$OUTPUT_DIR"

# Transcription factors to process
TRANSCRIPTION_FACTORS=("NANOG_ESC_1" "MYC_ESC_1" "OCT4_ESC_1" "YAP1_ESC_1" "NANOG_ESC_2" "MYC_ESC_2" "OCT4_ESC_2" "YAP1_ESC_2" "NANOG_EpiLC_1" "MYC_EpiLC_1" "OCT4_EpiLC_1" "YAP1_EpiLC_1" "NANOG_EpiLC_2" "MYC_EpiLC_2" "OCT4_EpiLC_2" "YAP1_EpiLC_2")

# Load samtools module
module load bio/SAMtools/1.19.2-GCC-13.2.0

# Log file for debugging
LOG_FILE="${OUTPUT_DIR}/merge_log.txt"
echo "Starting BAM file merge at $(date)" > "$LOG_FILE"

# Iterate over transcription factors
for TF in "${TRANSCRIPTION_FACTORS[@]}"; do
    # Create output file name for the current TF
    OUTPUT_FILE="${OUTPUT_DIR}/merged_sample_${TF}_mm39.bam"

    # Find all BAM files matching the transcription factor name
    bam_files_to_merge=("${BAM_DIR}"/R*C*_"${TF}"_*_bowtie2.mapped.bam)

    # Check if there are files to merge
    if [ ${#bam_files_to_merge[@]} -eq 0 ]; then
        echo "No BAM files found for transcription factor: $TF. Skipping." >> "$LOG_FILE"
        continue
    fi

    # Log the files to be merged
    echo "Merging ${#bam_files_to_merge[@]} files for $TF:" >> "$LOG_FILE"
    printf "%s\n" "${bam_files_to_merge[@]}" >> "$LOG_FILE"

    # Perform the merge
    samtools merge -@ 8 "$OUTPUT_FILE" "${bam_files_to_merge[@]}"

    # Check if the merge was successful
    if [ $? -eq 0 ]; then
        echo "Successfully merged BAM files for $TF into $OUTPUT_FILE" >> "$LOG_FILE"

        # Sort the merged BAM file
        SORTED_OUTPUT_FILE="${OUTPUT_DIR}/merged_sample_${TF}_mm39.sorted.bam"
        samtools sort -@ 8 -o "$SORTED_OUTPUT_FILE" "$OUTPUT_FILE"

        # Check if sorting was successful
        if [ $? -eq 0 ]; then
            echo "Successfully sorted merged BAM file for $TF into $SORTED_OUTPUT_FILE" >> "$LOG_FILE"
            # Remove the unsorted file to save space (optional)
            rm "$OUTPUT_FILE"
        else
            echo "Error sorting BAM file for $TF. Check samtools output." >> "$LOG_FILE"
        fi
    else
        echo "Error merging BAM files for $TF. Check samtools output." >> "$LOG_FILE"
    fi

done

echo "BAM file merging completed at $(date)" >> "$LOG_FILE"
```
## Merge single-cell BAM files into aggregated BAM file RHH mm10
```bash
nano aggregated_bam_mm10.sh

#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Directory containing BAM files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/bam_mm10"

# Output directory for merged BAM files
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm10"
mkdir -p "$OUTPUT_DIR"

# Transcription factors to process
TRANSCRIPTION_FACTORS=("NANOG_ESC_1" "MYC_ESC_1" "OCT4_ESC_1" "YAP1_ESC_1" "NANOG_ESC_2" "MYC_ESC_2" "OCT4_ESC_2" "YAP1_ESC_2" "NANOG_EpiLC_1" "MYC_EpiLC_1" "OCT4_EpiLC_1" "YAP1_EpiLC_1" "NANOG_EpiLC_2" "MYC_EpiLC_2" "OCT4_EpiLC_2" "YAP1_EpiLC_2")

# Load samtools module
module load bio/SAMtools/1.19.2-GCC-13.2.0

# Log file for debugging
LOG_FILE="${OUTPUT_DIR}/merge_log.txt"
echo "Starting BAM file merge at $(date)" > "$LOG_FILE"

# Iterate over transcription factors
for TF in "${TRANSCRIPTION_FACTORS[@]}"; do
    # Create output file name for the current TF
    OUTPUT_FILE="${OUTPUT_DIR}/merged_sample_${TF}_mm10.bam"

    # Find all BAM files matching the transcription factor name
    bam_files_to_merge=("${BAM_DIR}"/R*C*_"${TF}"_*_bowtie2.mapped.bam)

    # Check if there are files to merge
    if [ ${#bam_files_to_merge[@]} -eq 0 ]; then
        echo "No BAM files found for transcription factor: $TF. Skipping." >> "$LOG_FILE"
        continue
    fi

    # Log the files to be merged
    echo "Merging ${#bam_files_to_merge[@]} files for $TF:" >> "$LOG_FILE"
    printf "%s\n" "${bam_files_to_merge[@]}" >> "$LOG_FILE"

    # Perform the merge
    samtools merge -@ 8 "$OUTPUT_FILE" "${bam_files_to_merge[@]}"

    # Check if the merge was successful
    if [ $? -eq 0 ]; then
        echo "Successfully merged BAM files for $TF into $OUTPUT_FILE" >> "$LOG_FILE"

        # Sort the merged BAM file
        SORTED_OUTPUT_FILE="${OUTPUT_DIR}/merged_sample_${TF}_mm10.sorted.bam"
        samtools sort -@ 8 -o "$SORTED_OUTPUT_FILE" "$OUTPUT_FILE"

        # Check if sorting was successful
        if [ $? -eq 0 ]; then
            echo "Successfully sorted merged BAM file for $TF into $SORTED_OUTPUT_FILE" >> "$LOG_FILE"
            # Remove the unsorted file to save space (optional)
            rm "$OUTPUT_FILE"
        else
            echo "Error sorting BAM file for $TF. Check samtools output." >> "$LOG_FILE"
        fi
    else
        echo "Error merging BAM files for $TF. Check samtools output." >> "$LOG_FILE"
    fi

done

echo "BAM file merging completed at $(date)" >> "$LOG_FILE"
```
## Peak Calling via MACS2 from aggregated BAM files mm39
```bash
nano peak.calling_aggregated_bam_mm39.sh

#!/bin/bash
#SBATCH --job-name=PeakCalling
#SBATCH --output=peak_calling_%j.out
#SBATCH --error=peak_calling_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Directory containing sorted BAM files
INPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm39"

# Output directory for MACS2 peaks
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm39/peaks"
mkdir -p "$OUTPUT_DIR"

# Load required modules
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Perform peak calling for each sorted BAM file
for f in "$INPUT_DIR"/*sorted.bam; do
    BASENAME=$(basename "$f" .sorted.bam)
    macs2 callpeak -t "$f" -f BAMPE -g mm --keep-dup all -n "${OUTPUT_DIR}/${BASENAME}" \
    --nomodel --extsize 55 -B --SPMR
    if [ $? -eq 0 ]; then
        echo "MACS2 peak calling successful for $f"
    else
        echo "Error in MACS2 peak calling for $f" >&2
        exit 1
    fi
done

# Make 3-col bed file after peak calling

for f in "$OUTPUT_DIR"/*_peaks.narrowPeak; do
    BASENAME=$(basename "$f" _peaks.narrowPeak)
    awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - > "$OUTPUT_DIR/${BASENAME}_real.bed"
    if [ $? -eq 0 ]; then
        echo "3-column BED file created for $f"
    else
        echo "Error creating BED file for $f" >&2
        exit 1
    fi
done
```
## Peak Calling via MACS2 from aggregated BAM files mm10
```bash
nano peak.calling_aggregated_bam_mm10.sh

#!/bin/bash
#SBATCH --job-name=PeakCalling
#SBATCH --output=peak_calling_%j.out
#SBATCH --error=peak_calling_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Directory containing sorted BAM files
INPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm10"

# Output directory for MACS2 peaks
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/snDynaTag/fastq_20241203_GP_scDynaTag_mESC_EpiLC/outs/fastq_path/TRIMMED_DIR/alignment/aggregated_bam_mm10/peaks"
mkdir -p "$OUTPUT_DIR"

# Load required modules
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Perform peak calling for each sorted BAM file
for f in "$INPUT_DIR"/*sorted.bam; do
    BASENAME=$(basename "$f" .sorted.bam)
    macs2 callpeak -t "$f" -f BAMPE -g mm --keep-dup all -n "${OUTPUT_DIR}/${BASENAME}" \
    --nomodel --extsize 55 -B --SPMR
    if [ $? -eq 0 ]; then
        echo "MACS2 peak calling successful for $f"
    else
        echo "Error in MACS2 peak calling for $f" >&2
        exit 1
    fi
done

# Make 3-col bed file after peak calling

for f in "$OUTPUT_DIR"/*_peaks.narrowPeak; do
    BASENAME=$(basename "$f" _peaks.narrowPeak)
    awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - > "$OUTPUT_DIR/${BASENAME}_real.bed"
    if [ $? -eq 0 ]; then
        echo "3-column BED file created for $f"
    else
        echo "Error creating BED file for $f" >&2
        exit 1
    fi
done
```
