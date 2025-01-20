# Data Preprocessing
## Alignment mm39
```bash
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
```
## Alignment mm10
```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

module load bio/Bowtie2/2.5.1-GCC-12.3.0

# Reference genome index
ref='/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/ref/ref'

# Input and output directories
input_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/fastq'
output_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/sam'

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop over input files
for f1 in "$input_dir"/*_R1_001.fastq.gz; do
    # Construct the paired R2 file
    f2=${f1%_R1_001.fastq.gz}_R2_001.fastq.gz

    # Extract the base name (without path and extensions)
    base_name=$(basename "$f1" _R1_001.fastq.gz)

    # Construct the output SAM file path
    output_file="$output_dir/${base_name}_bowtie2.sam"

    # Run Bowtie2
    bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 \
        -x "$ref" -1 "$f1" -2 "$f2" -S "$output_file"
done
```

## Alignment mm10 GSE224292_mESC_CUTnRUN
```bash
#!/bin/bash

# SLURM parameters for each job
TIME="24:00:00"
MEM="16gb"
CPUS="8"
EMAIL="rhaensel@uni-koeln.de"

# Reference genome index
ref='/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/ref/ref'

# Input and output directories
input_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/fastq/GSE224292_mESC_CUTnRUN_fastq_clean'
output_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/sam/GSE224292_mESC_CUTnRUN_sam'

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop over input files
for f1 in "$input_dir"/*_1.fastq; do
    # Construct the paired R2 file
    f2=${f1%_1.fastq}_2.fastq

    # Extract the base name (without path and extensions)
    base_name=$(basename "$f1" _1.fastq)

    # Construct the output SAM file path
    output_file="$output_dir/${base_name}_bowtie2.sam"

    # Create a temporary SLURM script for each job
    job_script=$(mktemp)

    cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=bowtie2_${base_name}
#SBATCH --time=$TIME
#SBATCH --mem=$MEM
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$EMAIL

module load bio/Bowtie2/2.5.1-GCC-12.3.0

# Run Bowtie2
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p $CPUS -I 10 -X 700 \
    -x "$ref" -1 "$f1" -2 "$f2" -S "$output_file"
EOF

    # Submit the job
    sbatch "$job_script"

    # Optionally remove the temporary script (keep for debugging if needed)
    rm "$job_script"
done
```
## Alignment mm10 GSM4291125_mESC_ChIPseq
```bash
#!/bin/bash

# SLURM parameters for each job
TIME="24:00:00"
MEM="32gb"
CPUS="8"
EMAIL="rhaensel@uni-koeln.de"

# Reference genome index
ref='/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/ref/ref'

# Input and output directories
input_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/fastq/GSM4291125_mESC_ChIPseq_fastq_clean'
output_dir='/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/sam/GSM4291125_mESC_ChIPseq_sam'

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop over input files
for f1 in "$input_dir"/*_1.fastq; do
    # Construct the paired R2 file
    f2=${f1%_1.fastq}_2.fastq

    # Extract the base name (without path and extensions)
    base_name=$(basename "$f1" _1.fastq)

    # Construct the output SAM file path
    output_file="$output_dir/${base_name}_bowtie2.sam"

    # Create a temporary SLURM script for each job
    job_script=$(mktemp)

    cat > "$job_script" <<EOF
#!/bin/bash
#SBATCH --job-name=bowtie2_${base_name}
#SBATCH --time=$TIME
#SBATCH --mem=$MEM
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$EMAIL

module load bio/Bowtie2/2.5.1-GCC-12.3.0

# Run Bowtie2
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p $CPUS -I 10 -X 700 \
    -x "$ref" -1 "$f1" -2 "$f2" -S "$output_file"
EOF

    # Submit the job
    sbatch "$job_script"

    # Optionally remove the temporary script (keep for debugging if needed)
    rm "$job_script"
done
```
## Generate BAM files mm39
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 8GB --wrap "module load samtools/1.13 && samtools view -bS -F 0x04 "$f1" > /scratch/phunold/TFCT/ESC_d2EpiLC_ONSMY/alignment/bam/${f1%%}_bowtie2.mapped.bam"
done
```
## Generate BAM files mm10
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 32GB --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -bS -F 0x04 "$f1" > /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/${f1%%}_bowtie2.mapped.bam"
done
```
## Generate BAM files mm10 GSE224292_mESC_CUTnRUN
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 32GB --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -bS -F 0x04 "$f1" > /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSE224292_mESC_CUTnRUN_bam/${f1%%}_bowtie2.mapped.bam"
done
```
## Generate BAM files mm10 GSM4291125_mESC_ChIPseq
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 32GB --cpus-per-task 8 --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -bS -F 0x04 "$f1" > /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSM4291125_mESC_ChIPseq_bam/${f1%%}_bowtie2.mapped.bam"
done
```
## Sort and Remove duplicates files mm10 GSM4291125_mESC_ChIPseq
#!/bin/bash

# Directories
input_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSM4291125_mESC_ChIPseq_bam"
output_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSM4291125_mESC_ChIPseq_bam"

module load bio/SAMtools/1.19.2-GCC-13.2.0

for bam_file in "$input_dir"/*.bam; do
    base_name=$(basename "$bam_file" .bam)

    sbatch -J rmdup --mem=32GB --cpus-per-task=8 --wrap "
    module load bio/SAMtools/1.19.2-GCC-13.2.0 && \
    samtools sort -@ 8 -n \"$bam_file\" -o ${output_dir}/${base_name}_querysorted.bam && \
    samtools fixmate -@ 8 -m ${output_dir}/${base_name}_querysorted.bam ${output_dir}/${base_name}_fixmate.bam && \
    samtools sort -@ 8 ${output_dir}/${base_name}_fixmate.bam -o ${output_dir}/${base_name}_coordsorted.bam && \
    samtools markdup -@ 8 -r ${output_dir}/${base_name}_coordsorted.bam ${output_dir}/${base_name}_rmdup.bam && \
    rm ${output_dir}/${base_name}_querysorted.bam ${output_dir}/${base_name}_fixmate.bam ${output_dir}/${base_name}_coordsorted.bam"
done
```
## Downsampling mm39
```bash
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
```
## Downsampling mm10
```bash
#!/bin/bash

# Step 1: Count reads in BAM files
echo "Step 1: Counting reads in BAM files..."
num_bam_files=$(ls -1 *.bam | wc -l)
for f1 in *.bam; do
    if [[ ! -f "${f1%%.bam}.stat5" ]]; then
        sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=rhaensel@uni-koeln.de --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -c -F 260 $f1 > ${f1%%.bam}.stat5"
    fi
done

# Wait until all .stat5 files are generated
echo "Waiting for read count jobs to complete..."
while [ $(ls -1 *.stat5 2>/dev/null | wc -l) -lt $num_bam_files ]; do
    sleep 60
done
echo "Finished counting reads."

# Step 2: Combine results
echo "Step 2: Combining results..."
if [[ -f reads.txt ]]; then
    rm reads.txt  # Remove existing reads.txt to avoid appending duplicates
fi
cat *.stat5 > reads.txt
echo "Finished combining results."

# Step 3: Downsample BAM files
echo "Step 3: Downsampling BAM files..."
for file in *.bam; do
    # Determine desired reads based on file name
    if [[ "$file" == *"ATAC"* ]]; then
        desired_reads=15000000
    elif [[ "$file" == *"MYC"* || "$file" == *"YAP1"* ]]; then
        desired_reads=5000000
    elif [[ "$file" == *"OCT4"* ]]; then
        desired_reads=12000000
    elif [[ "$file" == *"SOX2"* ]]; then
        desired_reads=19000000
    elif [[ "$file" == *"H3K27me3"* || "$file" == *"H3K4me3"* || "$file" == *"CTCF"* || "$file" == *"NANOG"* ]]; then
        desired_reads=10000000
    else
        desired_reads=10000000  # Default desired reads
    fi

    # Extract base name (adjust according to your file naming convention)
    base_name="${file%_bowtie2.mapped.bam}"

    # Submit sbatch job
    sbatch --mem 8G --cpus-per-task 8 --wrap "
    module load bio/SAMtools/1.19.2-GCC-13.2.0
    total_reads=\$(samtools view -@ 8 -c -F 260 \"$file\")
    if [[ \$total_reads -ge $desired_reads ]]; then
        scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\")
        output_file=\"${base_name}_${desired_reads}_mm10_norm_clean.bam\"
        samtools view -@ 8 -bs \"\$scaling_factor\" \"$file\" > \"\$output_file\"
    else
        output_file=\"${base_name}_\${total_reads}_mm10_same_clean.bam\"
        cp \"$file\" \"\$output_file\"
    fi
    "
done

# Optionally, wait for all downsampling jobs to complete
echo "Waiting for downsampling jobs to complete..."
total_jobs=$(squeue -u $USER | grep "sbatch" | wc -l)
while [ $total_jobs -gt 0 ]; do
    sleep 60
    total_jobs=$(squeue -u $USER | grep "sbatch" | wc -l)
done
echo "Finished downsampling BAM files."

# remove redundancy in file name
rename '_L001_bowtie2.sam' '' *clean.bam
```
## Downsampling mm10 GSE224292_mESC_CUTnRUN
```bash
#!/bin/bash

# Step 1: Count reads in BAM files
echo "Step 1: Counting reads in BAM files..."
num_bam_files=$(ls -1 *.bam | wc -l)
for f1 in *.bam; do
    if [[ ! -f "${f1%%.bam}.stat5" ]]; then
        sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=rhaensel@uni-koeln.de --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -c -F 260 $f1 > ${f1%%.bam}.stat5"
    fi
done

# Wait until all .stat5 files are generated
echo "Waiting for read count jobs to complete..."
while [ $(ls -1 *.stat5 2>/dev/null | wc -l) -lt $num_bam_files ]; do
    sleep 60
done
echo "Finished counting reads."

# Step 2: Combine results
echo "Step 2: Combining results..."
if [[ -f reads.txt ]]; then
    rm reads.txt  # Remove existing reads.txt to avoid appending duplicates
fi
cat *.stat5 > reads.txt
echo "Finished combining results."

# Step 3: Downsample BAM files
echo "Step 3: Downsampling BAM files..."
for file in *.bam; do
    # Determine desired reads based on file name
    if [[ "$file" == *"SRR"* ]]; then
        desired_reads=10000000
   else
        desired_reads=10000000  # Default desired reads
    fi

    # Extract base name (adjust according to your file naming convention)
    base_name="${file%_bowtie2.mapped.bam}"

    # Submit sbatch job
    sbatch --mem 8G --cpus-per-task 8 --wrap "
    module load bio/SAMtools/1.19.2-GCC-13.2.0
    total_reads=\$(samtools view -@ 8 -c -F 260 \"$file\")
    if [[ \$total_reads -ge $desired_reads ]]; then
        scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\")
        output_file=\"${base_name}_${desired_reads}_mm10_norm_clean.bam\"
        samtools view -@ 8 -bs \"\$scaling_factor\" \"$file\" > \"\$output_file\"
    else
        output_file=\"${base_name}_\${total_reads}_mm10_same_clean.bam\"
        cp \"$file\" \"\$output_file\"
    fi
    "
done
```
## Downsampling mm10 GSM4291125_mESC_ChIPseq
```bash
#!/bin/bash

# Step 1: Count reads in BAM files
echo "Step 1: Counting reads in BAM files..."
num_bam_files=$(ls -1 *rmdup.bam | wc -l)
for f1 in *rmdup.bam; do
    if [[ ! -f "${f1%%.bam}.stat5" ]]; then
        sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=rhaensel@uni-koeln.de --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -c -F 260 $f1 > ${f1%%.bam}.stat5"
    fi
done

# Step 2: Combine results
echo "Step 2: Combining results..."
if [[ -f reads.txt ]]; then
    rm reads.txt  # Remove existing reads.txt to avoid appending duplicates
fi
cat *.stat5 > reads.txt
echo "Finished combining results."

# Step 3: Downsample BAM files
echo "Step 3: Downsampling BAM files..."
for file in *rmdup.bam; do
    # Determine desired reads based on file name
    if [[ "$file" == *"SRR"* ]]; then
        desired_reads=30000000
   else
        desired_reads=30000000  # Default desired reads
    fi

    # Extract base name (adjust according to your file naming convention)
    base_name="${file%_bowtie2.mapped_rmdup.bam}"

    # Submit sbatch job
    sbatch --mem 8G --cpus-per-task 8 --wrap "
    module load bio/SAMtools/1.19.2-GCC-13.2.0
    total_reads=\$(samtools view -@ 8 -c -F 260 \"$file\")
    if [[ \$total_reads -ge $desired_reads ]]; then
        scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\")
        output_file=\"${base_name}_${desired_reads}_mm10_norm_clean.bam\"
        samtools view -@ 8 -bs \"\$scaling_factor\" \"$file\" > \"\$output_file\"
    else
        output_file=\"${base_name}_\${total_reads}_mm10_same_clean.bam\"
        cp \"$file\" \"\$output_file\"
    fi
    "
done
# Step 4: Count reads in downsampled BAM files
echo "Step 4: Counting reads in downsampled BAM files..."
num_bam_files=$(ls -1 *_30000000_mm10_norm_clean.bam | wc -l)
for f1 in *_30000000_mm10_norm_clean.bam; do
    if [[ ! -f "${f1%%.bam}.stat5" ]]; then
        sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=rhaensel@uni-koeln.de --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools view -c -F 260 $f1 > ${f1%%.bam}.stat5"
    fi
done
# Step 5: Combine results
echo "Step 5: Combining results..."
if [[ -f reads_downsampled.txt ]]; then
    rm reads_downsampled.txt  # Remove existing reads.txt to avoid appending duplicates
fi
cat *_clean.stat5 > reads_downsampled.txt
echo "Finished combining results."
```
## Sort BAM files mm39
```bash
for f in *down_bowtie2.mapped.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load samtools/1.13 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
```
## Sort BAM files mm10
```bash
for f in *clean.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
```
## Sort BAM files mm10 GSE224292_mESC_CUTnRUN
```bash
for f in *clean.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
rename 'bowtie2.sam_' '' *_bowtie2.sam_*
rename 'SRR23310225' 'SRR23310225_SOX2_CUTnRUN_SL' *
rename 'SRR23310227' 'SRR23310227_SOX2_CUTnRUN_2i' *
rename 'SRR23310247' 'SRR23310247_OCT4_CUTnRUN_SL' *
rename 'SRR23310248' 'SRR23310248_OCT4_CUTnRUN_2i' *
rename 'SRR23310249' 'SRR23310249_NANOG_CUTnRUN_SL' *
rename 'SRR23310250' 'SRR23310250_NANOG_CUTnRUN_2i' *
rename 'SRR23310254' 'SRR23310254_IgG_CUTnRUN_SL' *
rename 'SRR23310255' 'SRR23310255_IgG_CUTnRUN_2i' *
```
## Sort and rename BAM files mm10 GSM4291125_mESC_ChIPseq
```bash
for f in *clean.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load bio/SAMtools/1.19.2-GCC-13.2.0 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
rename 'bowtie2.sam_' '' *_bowtie2.sam_*
rename 'SRR10992264' 'SRR10992264_WT_input_ChIPseq' *
rename 'SRR10992265' 'SRR10992265_NANOG_ChIPseq' *
rename 'SRR10992266' 'SRR10992266_OCT4_ChIPseq' *
rename 'SRR10992267' 'SRR10992267_SOX2_ChIPseq' *
rename 'SRR10992268' 'SRR10992268_YAP1_ChIPseq' *
```
## Call Peaks with MACS2 mm39
```bash
for f in *.sort.bam; do
sbatch -J MACS2 --mem 32GB --wrap "module load use.own && module load pypack/macs2 && macs2 callpeak -t $f -f BAMPE -g mm --keep-dup all -n /scratch/phunold/ESC/peaks/$f --nomodel --extsize 55 -B --SPMR"
done

module load bedtools/2.29.2
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
```
## Call Peaks with MACS2 mm10
```bash
for f in *.sort.bam; do
sbatch -J MACS2 --mem 32GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && macs2 callpeak -t $f -f BAMPE -g mm --keep-dup all -n /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/$f --nomodel --extsize 55 -B --SPMR"
done

for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
```
## Call Peaks with MACS2 mm10 GSE224292_mESC_CUTnRUN
```bash
#create script that considers matched IgG control as -c during peak calling 
nano run_macs2_callpeak.sh
./run_macs2_callpeak.sh

#!/bin/bash

# Define the mapping of treatment (-t) files to control (-c) files
declare -A controls
controls=(
    ["SRR23310225_SOX2_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"]="SRR23310254_IgG_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"
    ["SRR23310227_SOX2_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"]="SRR23310255_IgG_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"
    ["SRR23310247_OCT4_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"]="SRR23310254_IgG_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"
    ["SRR23310248_OCT4_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"]="SRR23310255_IgG_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"
    ["SRR23310249_NANOG_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"]="SRR23310254_IgG_CUTnRUN_SL_10000000_mm10_norm_clean.sort.bam"
    ["SRR23310250_NANOG_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"]="SRR23310255_IgG_CUTnRUN_2i_10000000_mm10_norm_clean.sort.bam"
)

# Directory to store peak results
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN"

# Loop over all treatment files and run MACS2 with the corresponding control
for t_file in *_mm10_norm_clean.sort.bam; do
    # Get the corresponding control file from the mapping
    c_file=${controls[$t_file]}
    
    # Check if the control file is found in the mapping
    if [ -z "$c_file" ]; then
        echo "No control file found for $t_file. Skipping..."
        continue
    fi

    # Construct the output prefix
    output_prefix="$peak_dir/$(basename "$t_file" .bam)"
    
    # Submit the SLURM job
    sbatch -J MACS2 --mem 32GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && macs2 callpeak -t $t_file -c $c_file -f BAMPE -g mm --keep-dup all -n $output_prefix --nomodel --extsize 55 -B --SPMR"
done

(abc-model-env) [rhaensel@ramses1 peaks_GSE224292_mESC_CUTnRUN]$ wc -l *narrowPeak
   1342 SRR23310225_SOX2_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
   1753 SRR23310227_SOX2_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak
     40 SRR23310247_OCT4_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
    299 SRR23310248_OCT4_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak
   8397 SRR23310249_NANOG_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
  19361 SRR23310250_NANOG_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak

# perform peak callign as exactly as done for DynaTag, don't consider IgG controls
nano run_macs2_callpeak_no_control.sh

#!/bin/bash

# Define the output directory
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN_no.control"

# Create the output directory if it doesn't exist
mkdir -p "$peak_dir"

# Loop over all treatment files and exclude IgG libraries
for t_file in *_mm10_norm_clean.sort.bam; do
    # Skip IgG libraries
    if [[ "$t_file" == *"IgG"* ]]; then
        echo "Skipping IgG library: $t_file"
        continue
    fi

    # Construct the output prefix
    output_prefix="$peak_dir/$(basename "$t_file" .bam)"
    
    # Submit the SLURM job
    sbatch -J MACS2 --mem 32GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && macs2 callpeak -t $t_file -f BAMPE -g mm --keep-dup all -n $output_prefix --nomodel --extsize 55 -B --SPMR"
done

(/home/rhaensel/intervene_env) [rhaensel@ramses1 peaks_GSE224292_mESC_CUTnRUN_no.control]$ wc -l *narrowPeak
   2318 SRR23310225_SOX2_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
   2479 SRR23310227_SOX2_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak
    216 SRR23310247_OCT4_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
    733 SRR23310248_OCT4_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak
   7242 SRR23310249_NANOG_CUTnRUN_SL_10000000_mm10_norm_clean.sort_peaks.narrowPeak
  15115 SRR23310250_NANOG_CUTnRUN_2i_10000000_mm10_norm_clean.sort_peaks.narrowPeak

/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN_no.control
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks_no.control.bed"
done
```
## Call Peaks with MACS2 mm10 GSM4291125_mESC_ChIPseq
```bash
#create script that considers matched Input control as -c during peak calling 
nano run_macs2_with_control.sh
chmod +x run_macs2_with_control.sh
./run_macs2_with_control.sh


#!/bin/bash

# Define the input control BAM file (used for all treatment files)
control_bam="SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean.sort.bam"

# Directory to store peak results
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks"

# Create the output directory if it does not exist
mkdir -p "$peak_dir"

# Loop over all treatment files and run MACS2 with the control
for t_file in *_30000000_mm10_norm_clean.sort.bam; do
    # Skip the control file
    if [ "$(basename "$t_file")" == "$control_bam" ]; then
        continue
    fi

    # Construct the output prefix
    output_prefix="$peak_dir/$(basename "$t_file" .bam)"

    # Submit the SLURM job
    sbatch -J MACS2 --mem 32GB --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    macs2 callpeak -t $t_file -c $control_bam -f BAMPE -g mm --keep-dup all -n $output_prefix --nomodel --extsize 55 -B --SPMR"
done

(/home/rhaensel/intervene_env) [rhaensel@ramses1 GSM4291125_mESC_ChIPseq_bam]$ wc -l ../../../peaks/GSM4291125_mESC_ChIPseq_peaks/*narrowPeak
   66848 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks/SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
   30512 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks/SRR10992266_OCT4_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
    5325 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks/SRR10992267_SOX2_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
    1686 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks/SRR10992268_YAP1_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak

# perform peak calling as exactly as done for DynaTag, don't consider Input
nano run_macs2_no_control.sh
chmod +x run_macs2_no_control.sh
./run_macs2_no_control.sh

#!/bin/bash

# Define the input control BAM file to exclude
control_bam="SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean.sort.bam"

# Directory to store peak results
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks_no_control"

# Create the output directory if it does not exist
mkdir -p "$peak_dir"

# Loop over all treatment files and run MACS2
for t_file in *_30000000_mm10_norm_clean.sort.bam; do
    # Skip the control file
    if [ "$(basename "$t_file")" == "$control_bam" ]; then
        continue
    fi

    # Construct the output prefix
    output_prefix="$peak_dir/$(basename "$t_file" .bam)"

    # Submit the SLURM job
    sbatch -J MACS2 --mem 32GB --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    macs2 callpeak -t $t_file -f BAMPE -g mm --keep-dup all -n $output_prefix --nomodel --extsize 55 -B --SPMR"
done

(/home/rhaensel/intervene_env) [rhaensel@ramses1 GSM4291125_mESC_ChIPseq_bam]$ wc -l ../../../peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/*narrowPeak
   73031 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
   37897 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/SRR10992266_OCT4_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
    8345 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/SRR10992267_SOX2_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak
    3573 ../../../peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/SRR10992268_YAP1_ChIPseq_30000000_mm10_norm_clean.sort_peaks.narrowPeak

cd /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
cd /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks_no_control
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks_no.control.bed"
done
```
## Calculate FRiP Scores mm39
```bash
#!/bin/bash

bam_dir="/scratch/phunold/ESC/bam"
peak_dir="/scratch/phunold/ESC/peaks"
output_file="FRiP_scores.txt"

module load samtools/1.13
module load bedtools/2.31.0

> "$output_file"

for bam_file in "$bam_dir"/*.sorted.bam; do
    sample_name=$(basename "$bam_file" .sorted.bam)
    peaks_bed="$peak_dir/${sample_name}_peaks.bed"
    
    if [ -f "$peaks_bed" ]; then
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $sample_name - FRiP Score: $FRiP" >> "$output_file"
    else
        echo "No matching peak file found for $sample_name"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Calculate FRiP Scores mm10 
```bash
nano FRiP_mm10_DynaTag.sh

#!/bin/bash

# Directories
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam"
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks"
output_file="FRiP_scores.txt"

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Initialize the output file
> "$output_file"

# Process each BAM file
for bam_file in "$bam_dir"/*_mm10_norm_clean.sort.bam; do
    # Extract the base name of the BAM file
    base_name=$(basename "$bam_file" _mm10_norm_clean.sort.bam)

    # Construct the corresponding peak file path
    peaks_bed="$peak_dir/${base_name}_mm10_norm_clean.sort.bam_peaks.narrowPeak_peaks.bed"

    # Check if the peak file exists
    if [ -f "$peaks_bed" ]; then
        echo "Processing BAM file: $bam_file with Peaks: $peaks_bed"

        # Calculate FRiP score
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $base_name - FRiP Score: $FRiP" >> "$output_file"
        
        # Clean up temporary files
        rm overlapping_reads.bam
    else
        echo "No matching peak file found for BAM file: $bam_file"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Calculate FRiP Scores mm10 GSE224292_mESC_CUTnRUN
```bash
nano FRiP_mm10_GSE224292_mESC_CUTnRUN.sh
chmod +x FRiP_mm10_GSE224292_mESC_CUTnRUN.sh

#!/bin/bash

# Directories
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSE224292_mESC_CUTnRUN_bam"
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN"
output_file="FRiP_scores.txt"

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Initialize the output file
> "$output_file"

# Process each BAM file
for bam_file in "$bam_dir"/*_10000000_mm10_norm_clean.sort.bam; do
    # Extract the base name of the BAM file
    base_name=$(basename "$bam_file" _10000000_mm10_norm_clean.sort.bam)

    # Construct the corresponding peak file path
    peaks_bed="$peak_dir/${base_name}_10000000_mm10_norm_clean.sort_peaks.narrowPeak_peaks.bed"

    # Check if the peak file exists
    if [ -f "$peaks_bed" ]; then
        echo "Processing BAM file: $bam_file with Peaks: $peaks_bed"

        # Calculate FRiP score
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $base_name - FRiP Score: $FRiP" >> "$output_file"
        
        # Clean up temporary files
        rm overlapping_reads.bam
    else
        echo "No matching peak file found for BAM file: $bam_file"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Calculate FRiP Scores mm10 GSE224292_mESC_CUTnRUN_no.control
```bash
nano FRiP_mm10_GSE224292_mESC_CUTnRUN_no_control.sh
#!/bin/bash

# Directories
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSE224292_mESC_CUTnRUN_bam"
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN_no.control"
output_file="FRiP_scores_no_control.txt"

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Initialize the output file
> "$output_file"

# Process each BAM file
for bam_file in "$bam_dir"/*_10000000_mm10_norm_clean.sort.bam; do
    # Extract the base name of the BAM file
    base_name=$(basename "$bam_file" _10000000_mm10_norm_clean.sort.bam)

    # Construct the corresponding peak file path
    peaks_bed="$peak_dir/${base_name}_10000000_mm10_norm_clean.sort_peaks.narrowPeak_peaks_no.control.bed"

    # Check if the peak file exists
    if [ -f "$peaks_bed" ]; then
        echo "Processing BAM file: $bam_file with Peaks: $peaks_bed"

        # Calculate FRiP score as a fraction
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $base_name - FRiP Score: $FRiP" >> "$output_file"

        # Clean up temporary files
        rm overlapping_reads.bam
    else
        echo "No matching peak file found for BAM file: $bam_file"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Calculate FRiP Scores mm10 GSM4291125_mESC_ChIPseq
```bash
nano FRiP_mm10_GSM4291125_mESC_ChIPseq.sh
chmod +x FRiP_mm10_GSM4291125_mESC_ChIPseq.sh
./FRiP_mm10_GSM4291125_mESC_ChIPseq.sh
#!/bin/bash

# Directories
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSM4291125_mESC_ChIPseq_bam"
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks"
output_file="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks/FRiP_scores_GSM4291125_mESC_ChIPseq.txt"

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Initialize the output file
> "$output_file"

# Process each BAM file
for bam_file in "$bam_dir"/*_30000000_mm10_norm_clean.sort.bam; do
    # Extract the base name of the BAM file
    base_name=$(basename "$bam_file" _30000000_mm10_norm_clean.sort.bam)

    # Construct the corresponding peak file path
    peaks_bed="$peak_dir/${base_name}_30000000_mm10_norm_clean.sort_peaks.narrowPeak_peaks.bed"

    # Check if the peak file exists
    if [ -f "$peaks_bed" ]; then
        echo "Processing BAM file: $bam_file with Peaks: $peaks_bed"

        # Calculate FRiP score as a fraction
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $base_name - FRiP Score: $FRiP" >> "$output_file"

        # Clean up temporary files
        rm overlapping_reads.bam
    else
        echo "No matching peak file found for BAM file: $bam_file"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Calculate FRiP Scores mm10 GSM4291125_mESC_ChIPseq_no.control
```bash
nano FRiP_mm10_GSM4291125_mESC_ChIPseq_no.control.sh
chmod +x FRiP_mm10_GSM4291125_mESC_ChIPseq_no.control.sh
./FRiP_mm10_GSM4291125_mESC_ChIPseq_no.control.sh
#!/bin/bash

# Directories
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/GSM4291125_mESC_ChIPseq_bam"
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks_no_control"
output_file="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks_no_control/FRiP_mm10_GSM4291125_mESC_ChIPseq_no.control.txt"

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Initialize the output file
> "$output_file"

# Process each BAM file
for bam_file in "$bam_dir"/*_30000000_mm10_norm_clean.sort.bam; do
    # Extract the base name of the BAM file
    base_name=$(basename "$bam_file" _30000000_mm10_norm_clean.sort.bam)

    # Construct the corresponding peak file path
    peaks_bed="$peak_dir/${base_name}_30000000_mm10_norm_clean.sort_peaks.narrowPeak_peaks_no.control.bed"

    # Check if the peak file exists
    if [ -f "$peaks_bed" ]; then
        echo "Processing BAM file: $bam_file with Peaks: $peaks_bed"

        # Calculate FRiP score as a fraction
        bedtools intersect -abam "$bam_file" -b "$peaks_bed" -wa -u > overlapping_reads.bam
        total_reads=$(samtools view -c -f 0x2 "$bam_file")
        reads_in_peaks=$(samtools view -c -f 0x2 overlapping_reads.bam)
        FRiP=$(bc -l <<< "$reads_in_peaks / $total_reads")
        echo "Sample: $base_name - FRiP Score: $FRiP" >> "$output_file"

        # Clean up temporary files
        rm overlapping_reads.bam
    else
        echo "No matching peak file found for BAM file: $bam_file"
    fi
done

echo "FRiP scores have been saved to $output_file"
```
## Filter Consensus Peaks mm39
```bash
module load bedtools/2.31.0

multiIntersectBed -i ESC-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G1_master_peaks.bed

multiIntersectBed -i ESC-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-S_master_peaks.bed

multiIntersectBed -i ESC-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G2_master_peaks.bed


multiIntersectBed -i ESC-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G1_master_peaks.bed

multiIntersectBed -i ESC-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-S_master_peaks.bed

multiIntersectBed -i ESC-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G2_master_peaks.bed


multiIntersectBed -i ESC-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G1_master_peaks.bed

multiIntersectBed -i ESC-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-S_master_peaks.bed

multiIntersectBed -i ESC-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G2_master_peaks.bed


multiIntersectBed -i ESC-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G1_master_peaks.bed

multiIntersectBed -i ESC-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-S_master_peaks.bed

multiIntersectBed -i ESC-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G2_master_peaks.bed


multiIntersectBed -i ESC-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G1_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G1_master_peaks.bed

multiIntersectBed -i ESC-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-S_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-S_master_peaks.bed

multiIntersectBed -i ESC-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G2_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G2_master_peaks.bed

multiIntersectBed -i ESC-ATAC*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-ATAC_master_peaks.bed
multiIntersectBed -i EpiLC-ATAC*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-ATAC_master_peaks.bed


cat ESC-OCT4-G1_master_peaks.bed EpiLC-OCT4-G1_master_peaks.bed | sortBed -i - > OCT4-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G1_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-S_master_peaks.bed EpiLC-OCT4-S_master_peaks.bed | sortBed -i - > OCT4-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-S_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-S_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-G2_master_peaks.bed EpiLC-OCT4-G2_master_peaks.bed | sortBed -i - > OCT4-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G2_all_peaks.over59nt.sorted.bed


cat ESC-SOX2-G1_master_peaks.bed EpiLC-SOX2-G1_master_peaks.bed | sortBed -i - > SOX2-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G1_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-S_master_peaks.bed EpiLC-SOX2-S_master_peaks.bed | sortBed -i - > SOX2-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-S_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-S_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-G2_master_peaks.bed EpiLC-SOX2-G2_master_peaks.bed | sortBed -i - > SOX2-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G2_all_peaks.over59nt.sorted.bed


cat ESC-NANOG-G1_master_peaks.bed EpiLC-NANOG-G1_master_peaks.bed | sortBed -i - > NANOG-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G1_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-S_master_peaks.bed EpiLC-NANOG-S_master_peaks.bed | sortBed -i - > NANOG-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-S_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-S_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-G2_master_peaks.bed EpiLC-NANOG-G2_master_peaks.bed | sortBed -i - > NANOG-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G2_all_peaks.over59nt.sorted.bed


cat ESC-YAP1-G1_master_peaks.bed EpiLC-YAP1-G1_master_peaks.bed | sortBed -i - > YAP1-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G1_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-S_master_peaks.bed EpiLC-YAP1-S_master_peaks.bed | sortBed -i - > YAP1-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-S_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-S_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-G2_master_peaks.bed EpiLC-YAP1-G2_master_peaks.bed | sortBed -i - > YAP1-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G2_all_peaks.over59nt.sorted.bed


cat ESC-MYC-G1_master_peaks.bed EpiLC-MYC-G1_master_peaks.bed | sortBed -i - > MYC-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G1_all_peaks.over59nt.sorted.bed

cat ESC-MYC-S_master_peaks.bed EpiLC-MYC-S_master_peaks.bed | sortBed -i - > MYC-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-S_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-S_all_peaks.over59nt.sorted.bed

cat ESC-MYC-G2_master_peaks.bed EpiLC-MYC-G2_master_peaks.bed | sortBed -i - > MYC-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G2_all_peaks.over59nt.sorted.bed


cat ESC-ATAC_master_peaks.bed EpiLC-ATAC_master_peaks.bed | sortBed -i - > ATAC_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' ATAC_all_peaks.bed | sortBed -i - | mergeBed -i - > ATAC_all_peaks.over59nt.sorted.bed
```
## Filter Consensus Peaks mm10 greater than 59bp peak size
```bash

nano filter_consensus_peaks_mm10.sh

multiIntersectBed -i ESC-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G1_master_peaks.bed

multiIntersectBed -i ESC-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-S_master_peaks.bed

multiIntersectBed -i ESC-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G2_master_peaks.bed

multiIntersectBed -i ESC-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G1_master_peaks.bed

multiIntersectBed -i ESC-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-S_master_peaks.bed

multiIntersectBed -i ESC-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G2_master_peaks.bed

multiIntersectBed -i ESC-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G1_master_peaks.bed

multiIntersectBed -i ESC-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-S_master_peaks.bed

multiIntersectBed -i ESC-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G2_master_peaks.bed

multiIntersectBed -i ESC-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G1_master_peaks.bed

multiIntersectBed -i ESC-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-S_master_peaks.bed

multiIntersectBed -i ESC-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G2_master_peaks.bed

multiIntersectBed -i ESC-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G1_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G1_master_peaks.bed

multiIntersectBed -i ESC-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-S_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-S_master_peaks.bed

multiIntersectBed -i ESC-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G2_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G2_master_peaks.bed

cat ESC-OCT4-G1_master_peaks.bed EpiLC-OCT4-G1_master_peaks.bed | sortBed -i - > OCT4-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G1_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-S_master_peaks.bed EpiLC-OCT4-S_master_peaks.bed | sortBed -i - > OCT4-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-S_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-S_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-G2_master_peaks.bed EpiLC-OCT4-G2_master_peaks.bed | sortBed -i - > OCT4-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G2_all_peaks.over59nt.sorted.bed


cat ESC-SOX2-G1_master_peaks.bed EpiLC-SOX2-G1_master_peaks.bed | sortBed -i - > SOX2-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G1_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-S_master_peaks.bed EpiLC-SOX2-S_master_peaks.bed | sortBed -i - > SOX2-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-S_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-S_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-G2_master_peaks.bed EpiLC-SOX2-G2_master_peaks.bed | sortBed -i - > SOX2-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G2_all_peaks.over59nt.sorted.bed


cat ESC-NANOG-G1_master_peaks.bed EpiLC-NANOG-G1_master_peaks.bed | sortBed -i - > NANOG-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G1_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-S_master_peaks.bed EpiLC-NANOG-S_master_peaks.bed | sortBed -i - > NANOG-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-S_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-S_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-G2_master_peaks.bed EpiLC-NANOG-G2_master_peaks.bed | sortBed -i - > NANOG-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G2_all_peaks.over59nt.sorted.bed


cat ESC-YAP1-G1_master_peaks.bed EpiLC-YAP1-G1_master_peaks.bed | sortBed -i - > YAP1-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G1_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-S_master_peaks.bed EpiLC-YAP1-S_master_peaks.bed | sortBed -i - > YAP1-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-S_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-S_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-G2_master_peaks.bed EpiLC-YAP1-G2_master_peaks.bed | sortBed -i - > YAP1-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G2_all_peaks.over59nt.sorted.bed


cat ESC-MYC-G1_master_peaks.bed EpiLC-MYC-G1_master_peaks.bed | sortBed -i - > MYC-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G1_all_peaks.over59nt.sorted.bed

cat ESC-MYC-S_master_peaks.bed EpiLC-MYC-S_master_peaks.bed | sortBed -i - > MYC-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-S_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-S_all_peaks.over59nt.sorted.bed

cat ESC-MYC-G2_master_peaks.bed EpiLC-MYC-G2_master_peaks.bed | sortBed -i - > MYC-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G2_all_peaks.over59nt.sorted.bed

rm *_all_peaks.bed
```
## Further Filter Consensus Peaks mm10 greater than 99bp peak size
```bash
nano filter_consensus_peaks_greater_than_99bp_mm10.sh

cat ESC-OCT4-G1_master_peaks.bed EpiLC-OCT4-G1_master_peaks.bed | sortBed -i - > OCT4-G1_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' OCT4-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G1_all_peaks.over99nt.sorted.bed

cat ESC-OCT4-S_master_peaks.bed EpiLC-OCT4-S_master_peaks.bed | sortBed -i - > OCT4-S_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' OCT4-S_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-S_all_peaks.over99nt.sorted.bed

cat ESC-OCT4-G2_master_peaks.bed EpiLC-OCT4-G2_master_peaks.bed | sortBed -i - > OCT4-G2_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' OCT4-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G2_all_peaks.over99nt.sorted.bed


cat ESC-SOX2-G1_master_peaks.bed EpiLC-SOX2-G1_master_peaks.bed | sortBed -i - > SOX2-G1_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' SOX2-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G1_all_peaks.over99nt.sorted.bed

cat ESC-SOX2-S_master_peaks.bed EpiLC-SOX2-S_master_peaks.bed | sortBed -i - > SOX2-S_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' SOX2-S_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-S_all_peaks.over99nt.sorted.bed

cat ESC-SOX2-G2_master_peaks.bed EpiLC-SOX2-G2_master_peaks.bed | sortBed -i - > SOX2-G2_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' SOX2-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G2_all_peaks.over99nt.sorted.bed


cat ESC-NANOG-G1_master_peaks.bed EpiLC-NANOG-G1_master_peaks.bed | sortBed -i - > NANOG-G1_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' NANOG-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G1_all_peaks.over99nt.sorted.bed

cat ESC-NANOG-S_master_peaks.bed EpiLC-NANOG-S_master_peaks.bed | sortBed -i - > NANOG-S_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' NANOG-S_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-S_all_peaks.over99nt.sorted.bed

cat ESC-NANOG-G2_master_peaks.bed EpiLC-NANOG-G2_master_peaks.bed | sortBed -i - > NANOG-G2_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' NANOG-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G2_all_peaks.over99nt.sorted.bed


cat ESC-YAP1-G1_master_peaks.bed EpiLC-YAP1-G1_master_peaks.bed | sortBed -i - > YAP1-G1_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' YAP1-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G1_all_peaks.over99nt.sorted.bed

cat ESC-YAP1-S_master_peaks.bed EpiLC-YAP1-S_master_peaks.bed | sortBed -i - > YAP1-S_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' YAP1-S_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-S_all_peaks.over99nt.sorted.bed

cat ESC-YAP1-G2_master_peaks.bed EpiLC-YAP1-G2_master_peaks.bed | sortBed -i - > YAP1-G2_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' YAP1-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G2_all_peaks.over99nt.sorted.bed


cat ESC-MYC-G1_master_peaks.bed EpiLC-MYC-G1_master_peaks.bed | sortBed -i - > MYC-G1_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' MYC-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G1_all_peaks.over99nt.sorted.bed

cat ESC-MYC-S_master_peaks.bed EpiLC-MYC-S_master_peaks.bed | sortBed -i - > MYC-S_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' MYC-S_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-S_all_peaks.over99nt.sorted.bed

cat ESC-MYC-G2_master_peaks.bed EpiLC-MYC-G2_master_peaks.bed | sortBed -i - > MYC-G2_all_peaks.bed
awk '{if (($3-$2) >= 99) print $0}' MYC-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G2_all_peaks.over99nt.sorted.bed

rm *_all_peaks.bed
```
## Further Filter Consensus Peaks mm10 greater than 149bp peak size
```bash
nano filter_consensus_peaks_greater_than_149bp_mm10.sh

cat ESC-OCT4-G1_master_peaks.bed EpiLC-OCT4-G1_master_peaks.bed | sortBed -i - > OCT4-G1_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' OCT4-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G1_all_peaks.over149nt.sorted.bed

cat ESC-OCT4-S_master_peaks.bed EpiLC-OCT4-S_master_peaks.bed | sortBed -i - > OCT4-S_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' OCT4-S_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-S_all_peaks.over149nt.sorted.bed

cat ESC-OCT4-G2_master_peaks.bed EpiLC-OCT4-G2_master_peaks.bed | sortBed -i - > OCT4-G2_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' OCT4-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G2_all_peaks.over149nt.sorted.bed


cat ESC-SOX2-G1_master_peaks.bed EpiLC-SOX2-G1_master_peaks.bed | sortBed -i - > SOX2-G1_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' SOX2-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G1_all_peaks.over149nt.sorted.bed

cat ESC-SOX2-S_master_peaks.bed EpiLC-SOX2-S_master_peaks.bed | sortBed -i - > SOX2-S_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' SOX2-S_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-S_all_peaks.over149nt.sorted.bed

cat ESC-SOX2-G2_master_peaks.bed EpiLC-SOX2-G2_master_peaks.bed | sortBed -i - > SOX2-G2_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' SOX2-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G2_all_peaks.over149nt.sorted.bed


cat ESC-NANOG-G1_master_peaks.bed EpiLC-NANOG-G1_master_peaks.bed | sortBed -i - > NANOG-G1_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' NANOG-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G1_all_peaks.over149nt.sorted.bed

cat ESC-NANOG-S_master_peaks.bed EpiLC-NANOG-S_master_peaks.bed | sortBed -i - > NANOG-S_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' NANOG-S_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-S_all_peaks.over149nt.sorted.bed

cat ESC-NANOG-G2_master_peaks.bed EpiLC-NANOG-G2_master_peaks.bed | sortBed -i - > NANOG-G2_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' NANOG-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G2_all_peaks.over149nt.sorted.bed


cat ESC-YAP1-G1_master_peaks.bed EpiLC-YAP1-G1_master_peaks.bed | sortBed -i - > YAP1-G1_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' YAP1-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G1_all_peaks.over149nt.sorted.bed

cat ESC-YAP1-S_master_peaks.bed EpiLC-YAP1-S_master_peaks.bed | sortBed -i - > YAP1-S_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' YAP1-S_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-S_all_peaks.over149nt.sorted.bed

cat ESC-YAP1-G2_master_peaks.bed EpiLC-YAP1-G2_master_peaks.bed | sortBed -i - > YAP1-G2_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' YAP1-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G2_all_peaks.over149nt.sorted.bed


cat ESC-MYC-G1_master_peaks.bed EpiLC-MYC-G1_master_peaks.bed | sortBed -i - > MYC-G1_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' MYC-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G1_all_peaks.over149nt.sorted.bed

cat ESC-MYC-S_master_peaks.bed EpiLC-MYC-S_master_peaks.bed | sortBed -i - > MYC-S_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' MYC-S_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-S_all_peaks.over149nt.sorted.bed

cat ESC-MYC-G2_master_peaks.bed EpiLC-MYC-G2_master_peaks.bed | sortBed -i - > MYC-G2_all_peaks.bed
awk '{if (($3-$2) >= 149) print $0}' MYC-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G2_all_peaks.over149nt.sorted.bed

rm *_all_peaks.bed
```
## Generate Count Matrices mm39
```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/phunold/ESC_EpiLC/bam/peaks/consensus/master_peaks"
bam_dir="/scratch/phunold/ESC_EpiLC/bam"
bedgraph_path="/scratch/phunold/ESC_EpiLC/bedgraph"

module load bedtools/2.29.2

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over59nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*.sorted.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%.sorted.bam}).bedgraph"
                    awk -v OFS='\t' '{print "chr"$1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%.sorted.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%.sorted.bam})_counts.txt"
                done
            fi
        done
    done
done
```
## Generate Count Matrices mm10 peak size over 59bp
```bash
nano Generate_Count_Matrices_mm10.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_mm10"

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over59nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_norm_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done

nano Generate_Count_Matrices_mm10_not.norm.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_mm10"
for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over59nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_same_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done
```
## Generate Count Matrices mm10 peak size over 99bp
```bash
nano Generate_Count_Matrices_over99bp_mm10.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_over99bp_mm10"

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over99nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_norm_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done

nano Generate_Count_Matrices_over99bp_mm10_not.norm.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_over99bp_mm10"
for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over99nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_same_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done
```
## Generate Count Matrices mm10 peak size over 149bp
```bash
nano Generate_Count_Matrices_over149bp_mm10.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_over149bp_mm10"

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over149nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_norm_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done

nano Generate_Count_Matrices_over149bp_mm10_not.norm.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ESC_EpiLC_bulk_DynaTag_CUTnTag_ATAC"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_over149bp_mm10"
for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over149nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_same_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_same_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done
```
## Generate BigWig Files for optical inspection mm39
```bash

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-ATAC-*_down_bowtie2.mapped.sort.bam | samtools sort -o EpiLC-ATAC.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-ATAC*_down_bowtie2.mapped.sort.bam | samtools sort -o ESC-ATAC.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-CTCF-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-CTCF-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-CTCF-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-CTCF-W300.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K27me3-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K27me3-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K27me3-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K27me3-W300.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K4me3-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K4me3-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K4me3-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K4me3-W300.merged.sorted.bam -"

for f1 in *merged.sorted.bam; do
sbatch -J BAMindex --mem 8GB --wrap "module load samtools/1.13 && samtools index $f1"
done

for f1 in *merged.sorted.bam; do
sbatch --mem 16G -J BAM2BW --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && bamCoverage -p 8 -b $f1 -bs 10 --skipNAs --centerReads --normalizeUsing CPM -of bigwig -o ./${f1%%.sorted.bam}_cpm.bw "
done
```
## Generate BigWig Files for optical inspection mm10
```bash

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-2d-MYC-G1-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-2d-MYC-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-2d-MYC-S-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-2d-MYC-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-2d-MYC-G2-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-2d-MYC-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-MYC-G1-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-MYC-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-MYC-S-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-MYC-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-MYC-G2-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-MYC-G2.mm10.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-OCT4-G1-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-OCT4-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-OCT4-S-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-OCT4-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-OCT4-G2-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-OCT4-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-OCT4-G1-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-OCT4-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-OCT4-S-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-OCT4-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-OCT4-G2-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-OCT4-G2.mm10.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-SOX2-G1-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-SOX2-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-SOX2-S-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-SOX2-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-SOX2-G2-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-SOX2-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-SOX2-G1-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-SOX2-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-SOX2-S-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-SOX2-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-SOX2-G2-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-SOX2-G2.mm10.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-NANOG-G1-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-NANOG-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-NANOG-S-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-NANOG-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-NANOG-G2-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-NANOG-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-NANOG-G1-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-NANOG-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-NANOG-S-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-NANOG-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-NANOG-G2-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-NANOG-G2.mm10.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-YAP1-G1-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-YAP1-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-YAP1-S-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-YAP1-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - EpiLC-d2-YAP1-G2-*_mm10_norm_clean.sort.bam | samtools sort -o EpiLC-d2-YAP1-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-YAP1-G1-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-YAP1-G1.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-YAP1-S-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-YAP1-S.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-YAP1-G2-*_mm10_norm_clean.sort.bam | samtools sort -o ESC-YAP1-G2.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-CTCF-ICS*__mm10_norm_clean.sort.bam | samtools sort -o ESC-CTCF-ICS.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-CTCF-W300*__mm10_norm_clean.sort.bam | samtools sort -o ESC-CTCF-W300.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-H3K27me3-ICS*__mm10_norm_clean.sort.bam | samtools sort -o ESC-H3K27me3-ICS.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-H3K27me3-W300*__mm10_norm_clean.sort.bam | samtools sort -o ESC-H3K27me3-W300.mm10.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-H3K4me3-ICS*__mm10_norm_clean.sort.bam | samtools sort -o ESC-H3K4me3-ICS.mm10.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools merge - ESC-H3K4me3-W300*__mm10_norm_clean.sort.bam | samtools sort -o ESC-H3K4me3-W300.mm10.merged.sorted.bam -"

for f1 in *merged.sorted.bam; do
sbatch -J BAMindex --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools index $f1"
done

for f1 in *merged.sorted.bam; do
sbatch --mem 16G -J BAM2BW --cpus-per-task 8 --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bamCoverage -p 8 -b $f1 -bs 10 --skipNAs --centerReads --normalizeUsing CPM -of bigwig -o /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/${f1%%.sorted.bam}_cpm.bw"
done
```
## Generate BigWig Files for optical inspection mm10 GSE224292_mESC_CUTnRUN --> not performed since spike-in normalised bigwigs are available under GSE224292 by authors and used for comparison with DynaTag mm10 and GSM4291125_mESC_ChIPseq mm10

## Generate BigWig Files for optical inspection mm10 GSM4291125_mESC_ChIPseq_bam
```bash

for f1 in *_mm10_norm_clean.sort.bam; do
sbatch -J BAMindex --mem 8GB --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && samtools index $f1"
done

for f1 in *_mm10_norm_clean.sort.bam; do
sbatch --mem 16G -J BAM2BW --cpus-per-task 8 --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bamCoverage -p 8 -b $f1 -bs 10 --skipNAs --centerReads --normalizeUsing CPM -of bigwig -o /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/${f1%%.sort.bam}_cpm.bw"
done
```
## plotCorrelation mm39
```bash

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-CTCF-W300.merged.sorted.bam ESC-CTCF-ICS.merged.sorted.bam -p 8 -o CTCF_W300vsICS.npz "
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-H3K27me3-W300.merged.sorted.bam ESC-H3K27me3-ICS.merged.sorted.bam -p 8 -o H3K27me3_W300vsICS.npz "
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-H3K4me3-W300.merged.sorted.bam ESC-H3K4me3-ICS.merged.sorted.bam -p 8 -o H3K4me3_W300vsICS.npz "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in CTCF_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o CTCF_W300vsICS_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in H3K27me3_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o H3K27me3_W300vsICS_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in H3K4me3_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o H3K4me3_W300vsICS_plotCorrelation.pdf "


sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-G1-1_S1_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-G1-2_S4_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-S-1_S2_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-S-2_S5_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-G2-1_S3_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-G2-2_S6_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-G1-1_S7_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-G1-2_S10_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-S-1_S8_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-S-2_S11_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-G2-1_S9_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-G2-2_S12_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-G1-1_S13_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-G1-2_S16_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-S-1_S14_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-S-2_S17_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-G2-1_S15_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-G2-2_S18_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-G1-1_S1_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-G1-2_S4_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-S-1_S2_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-S-2_S5_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-G2-1_S3_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-G2-2_S6_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-G1-1_S25_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-G1-2_S28_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-S-1_S26_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-S-2_S29_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-G2-1_S27_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-G2-2_S30_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_G2.npz"

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_G2_plotCorrelation.pdf "
```
## Inter Size Plot mm39
```bash
for f in *ICS*merged.sorted.bam; do
 sbatch --mem 8G --time=2:00:00 -J size --wrap "module load openjdk/1.8.0_60 && java -Xmx7g -jar /home/phunold/privatemodules/pypack/picard.jar CollectInsertSizeMetrics \
      I=$f \
      O=./${f%%}_size_metrics.txt \
      H=./${f%%}_size_histogram.pdf"
done

for f in *W300*merged.sorted.bam; do
 sbatch --mem 8G --time=2:00:00 -J size --wrap "module load openjdk/1.8.0_60 && java -Xmx7g -jar /home/phunold/privatemodules/pypack/picard.jar CollectInsertSizeMetrics \
      I=$f \
      O=./${f%%}_size_metrics.txt \
      H=./${f%%}_size_histogram.pdf"
done
```
# Differential Binding Analysis mm39
## Load Packages
```R
library(edgeR)
library(ggplot2)

```
## Load Data and Assemble Count Matrices
```R
rm(list = ls())

# Define the file path
file_path <- "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/count_matrices/"

# Define the TF to analyze
TF <- "OCT4"  # Change this to any other TF you want to analyze

# Read count data from files
ESC_G1_1 <- read.table(paste0(file_path, "ESC-", TF, "-G1-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_G1_2 <- read.table(paste0(file_path, "ESC-", TF, "-G1-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G1_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G1-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G1_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G1-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

ESC_S_1 <- read.table(paste0(file_path, "ESC-", TF, "-S-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_S_2 <- read.table(paste0(file_path, "ESC-", TF, "-S-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_S_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-S-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_S_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-S-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

ESC_G2_1 <- read.table(paste0(file_path, "ESC-", TF, "-G2-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_G2_2 <- read.table(paste0(file_path, "ESC-", TF, "-G2-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G2_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G2-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G2_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G2-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

# Check if row names are identical in all matrices for G1 phase
identical_row_names_G1 <- all(
  identical(rownames(ESC_G1_1), rownames(ESC_G1_2)),
  identical(rownames(ESC_G1_1), rownames(EpiLC_G1_1)),
  identical(rownames(ESC_G1_1), rownames(EpiLC_G1_2))
)

# Check if row names are identical in all matrices for S phase
identical_row_names_S <- all(
  identical(rownames(ESC_S_1), rownames(ESC_S_2)),
  identical(rownames(ESC_S_1), rownames(EpiLC_S_1)),
  identical(rownames(ESC_S_1), rownames(EpiLC_S_2))
)

# Check if row names are identical in all matrices for G2 phase
identical_row_names_G2 <- all(
  identical(rownames(ESC_G2_1), rownames(ESC_G2_2)),
  identical(rownames(ESC_G2_1), rownames(EpiLC_G2_1)),
  identical(rownames(ESC_G2_1), rownames(EpiLC_G2_2))
)


ESC_G1_1 <- ESC_G1_1
ESC_G1_2 <- ESC_G1_2[, 1]
EpiLC_G1_1 <- EpiLC_G1_1[, 1]
EpiLC_G1_2 <- EpiLC_G1_2[, 1]

ESC_S_1 <- ESC_S_1
ESC_S_2 <- ESC_S_2[, 1]
EpiLC_S_1 <- EpiLC_S_1[, 1]
EpiLC_S_2 <- EpiLC_S_2[, 1]

ESC_G2_1 <- ESC_G2_1
ESC_G2_2 <- ESC_G2_2[, 1]
EpiLC_G2_1 <- EpiLC_G2_1[, 1]
EpiLC_G2_2 <- EpiLC_G2_2[, 1]


# Merging Data
G1_merged <- cbind(ESC_G1_1, ESC_G1_2, EpiLC_G1_1, EpiLC_G1_2)
colnames(G1_merged)[1] <- "ESC_G1_1"

S_merged <- cbind(ESC_S_1, ESC_S_2, EpiLC_S_1, EpiLC_S_2)
colnames(S_merged)[1] <- "ESC_S_1"

G2_merged <- cbind(ESC_G2_1, ESC_G2_2, EpiLC_G2_1, EpiLC_G2_2)
colnames(G2_merged)[1] <- "ESC_G2_1"
```
## Setup edgeR Parametres
```R
# Set the reference group manually (change "PBS" to the group you want as the reference)
reference_group <- "ESC"

# Create a factor variable for group and specify the reference level
group <- factor(c("ESC", "ESC", "EpiLC", "EpiLC"), levels = c("ESC", "EpiLC"))

# Design Matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Create contrast with ESC as the reference group
contrast <- makeContrasts(reference_vs_other = paste0("EpiLC - ", reference_group), levels = design)
```
## Differential Binding Analysis per Cell Cycle Phase
```R
# Data Preprocessing
dge_G1 <- DGEList(counts = G1_merged)
dge_G1$samples$group <- group  # Adding group information to DGEList object
dge_G1 <- calcNormFactors(dge_G1, method = "TMM")

# Plot MDS
plotMDS(dge_G1, main="", col=as.numeric(group)) # Coloring by group

dge_G1 <- estimateGLMCommonDisp(dge_G1, design)
dge_G1 <- estimateGLMTagwiseDisp(dge_G1, design)

# Estimate Dispersion
dge_G1 <- estimateDisp(dge_G1)

# Fit GLM
fit_G1 <- glmFit(dge_G1, design)

# Likelihood Ratio Test (LRT)
lrt_G1 <- glmLRT(fit_G1, contrast = contrast)

# Extract and adjust results
results_G1 <- topTags(lrt_G1, n = Inf)
results_G1 <- results_G1$table
results_G1$FDR <- p.adjust(results_G1$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_G1)

# Volcano Plot
res_G1 <- na.omit(results_G1)
sig_G1 <- sum(res_G1$FDR < 0.05)
UP_G1 <- sum(res_G1$logFC > 0.5 & res_G1$FDR < 0.05)
DOWN_G1 <- sum(res_G1$logFC < -0.5 & res_G1$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_G1, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_G1, ")"))
}
ggplot(data = res_G1, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_G1$logFC)), max(abs(res_G1$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "G0/G1 Phase"))




# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_G1 <- DOWN_G1

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_G1 <- UP_G1

# Create a data frame for the counts
count_data_G1 <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_G1, gained_regions_G1)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_G1 <- count_data_G1 %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_G1 <- ggplot(count_data_G1, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in G0/G1 Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_G1)

ggsave(filename = paste0(file_path, TF, "_Regions_G1.pdf"), plot = plot_G1)



###### S #####

# Data Preprocessing
dge_S <- DGEList(counts = S_merged)
dge_S$samples$group <- group  # Adding group information to DGEList object
dge_S <- calcNormFactors(dge_S, method = "TMM")

# Plot MDS
plotMDS(dge_S, main="", col=as.numeric(group)) # Coloring by group

dge_S <- estimateGLMCommonDisp(dge_S, design)
dge_S <- estimateGLMTagwiseDisp(dge_S, design)

# Estimate Dispersion
dge_S <- estimateDisp(dge_S)

# Fit GLM
fit_S <- glmFit(dge_S, design)

# Likelihood Ratio Test (LRT)
lrt_S <- glmLRT(fit_S, contrast = contrast)

# Extract and adjust results
results_S <- topTags(lrt_S, n = Inf)
results_S <- results_S$table
results_S$FDR <- p.adjust(results_S$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_S)

# Volcano Plot
res_S <- na.omit(results_S)
sig_S <- sum(res_S$FDR < 0.05)
UP_S <- sum(res_S$logFC > 0.5 & res_S$FDR < 0.05)
DOWN_S <- sum(res_S$logFC < -0.5 & res_S$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_S, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_S, ")"))
}
ggplot(data = res_S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_S$logFC)), max(abs(res_S$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "S Phase"))
###

# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_S <- DOWN_S

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_S <- UP_S

# Create a data frame for the counts
count_data_S <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_S, gained_regions_S)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_S <- count_data_S %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_S <- ggplot(count_data_S, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in S Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_S)

ggsave(filename = paste0(file_path, TF, "_Regions_S.pdf"), plot = plot_S)



###### G2 #####

# Data Preprocessing
dge_G2 <- DGEList(counts = G2_merged)
dge_G2$samples$group <- group  # Adding group information to DGEList object
dge_G2 <- calcNormFactors(dge_G2, method = "TMM")

# Plot MDS
plotMDS(dge_G2, main="", col=as.numeric(group)) # Coloring by group

dge_G2 <- estimateGLMCommonDisp(dge_G2, design)
dge_G2 <- estimateGLMTagwiseDisp(dge_G2, design)

# Estimate Dispersion
dge_G2 <- estimateDisp(dge_G2)

# Fit GLM
fit_G2 <- glmFit(dge_G2, design)

# Likelihood Ratio Test (LRT)
lrt_G2 <- glmLRT(fit_G2, contrast = contrast)

# Extract and adjust results
results_G2 <- topTags(lrt_G2, n = Inf)
results_G2 <- results_G2$table
results_G2$FDR <- p.adjust(results_G2$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_G2)

#Volcano Plot
res_G2 <- na.omit(results_G2)
sig_G2 <- sum(res_G2$FDR < 0.05)
UP_G2 <- sum(res_G2$logFC > 0.5 & res_G2$FDR < 0.05)
DOWN_G2 <- sum(res_G2$logFC < -0.5 & res_G2$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_G2, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_G2, ")"))
}
ggplot(data = res_G2, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_G2$logFC)), max(abs(res_G2$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "G2/M Phase"))
###

# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_G2 <- DOWN_G2

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_G2 <- UP_G2

# Create a data frame for the counts
count_data_G2 <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_G2, gained_regions_G2)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_G2 <- count_data_G2 %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_G2 <- ggplot(count_data_G2, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in G2/M Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_G2)

ggsave(filename = paste0(file_path, TF, "_Regions_G2.pdf"), plot = plot_G2)




# Combine count data for all phases
combined_count_data <- rbind(count_data_G1, count_data_S, count_data_G2)

# Add a new column for phase
combined_count_data$Phase <- c(rep("G0/G1 Phase", nrow(count_data_G1)),
                               rep("S Phase", nrow(count_data_S)),
                               rep("G2/M Phase", nrow(count_data_G2)))

# Reorder the factor levels of Phase
combined_count_data$Phase <- factor(combined_count_data$Phase, levels = c("G0/G1 Phase", "S Phase", "G2/M Phase"))

# Reorder the factor levels of Category within each phase
combined_count_data <- combined_count_data %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

# Create the side-by-side bar chart with facets
plot <- ggplot(combined_count_data, aes(x = Phase, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of Differential", TF, "Binding"), x = "Cell Cycle Phase", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))

plot

# Save the plot as a PDF
ggsave(filename = paste0(file_path, TF, "_Regions_Cell_Cycle_Phases.pdf"), plot = plot)
```

# RHH reproduces Pascal's Differential Binding Analysis mm39
## Load Packages
```R
# RHH reproduces Pascal's Differential Binding Analysis mm39
## Load Packages
library(edgeR)
library(ggplot2)
library(dplyr)

## Load Data and Assemble Count Matrices
rm(list = ls())

# Define the file path
file_path <- "/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices/"

# Define the TF to analyze
TF <- "NANOG"  # Change this to any other TF you want to analyze

# Read count data from files for all phases
ESC_G1_1 <- read.table(paste0(file_path, "ESC-", TF, "-G1-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_G1_2 <- read.table(paste0(file_path, "ESC-", TF, "-G1-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G1_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G1-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G1_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G1-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

ESC_S_1 <- read.table(paste0(file_path, "ESC-", TF, "-S-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_S_2 <- read.table(paste0(file_path, "ESC-", TF, "-S-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_S_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-S-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_S_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-S-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

ESC_G2_1 <- read.table(paste0(file_path, "ESC-", TF, "-G2-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
ESC_G2_2 <- read.table(paste0(file_path, "ESC-", TF, "-G2-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G2_1 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G2-1_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)
EpiLC_G2_2 <- read.table(paste0(file_path, "EpiLC-d2-", TF, "-G2-2_counts.txt"), header=FALSE, row.names=1, check.names=FALSE)

# Prepare merged data
merge_phase_data <- function(esc1, esc2, epilc1, epilc2) {
  esc2 <- esc2[, 1]
  epilc1 <- epilc1[, 1]
  epilc2 <- epilc2[, 1]
  merged <- cbind(esc1, esc2, epilc1, epilc2)
  colnames(merged)[1] <- "ESC_1"
  return(merged)
}

G1_merged <- merge_phase_data(ESC_G1_1, ESC_G1_2, EpiLC_G1_1, EpiLC_G1_2)
S_merged <- merge_phase_data(ESC_S_1, ESC_S_2, EpiLC_S_1, EpiLC_S_2)
G2_merged <- merge_phase_data(ESC_G2_1, ESC_G2_2, EpiLC_G2_1, EpiLC_G2_2)

## Setup edgeR Parameters
analyze_phase <- function(merged_data, phase, file_path, TF) {
  # Set up design and group
  group <- factor(c("ESC", "ESC", "EpiLC", "EpiLC"), levels = c("ESC", "EpiLC"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrast <- makeContrasts(reference_vs_other = "EpiLC - ESC", levels = design)
  
  # edgeR Differential Binding Analysis
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Plot MDS
  pdf(file = paste0(file_path, TF, "_MDS_Plot_", phase, ".pdf"))
  plotMDS(dge, main = paste(TF, "MDS Plot -", phase), col = as.numeric(group))
  dev.off()
  
  # Fit the model
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract results
  results <- topTags(lrt, n = Inf)$table
  results$FDR <- p.adjust(results$PValue, method = "BH")
  
  # Volcano Plot with Correct Labels
  sig_regions <- results %>%
    mutate(
      Category = case_when(
        logFC < -0.5 & FDR < 0.05 ~ "Downregulated",
        logFC > 0.5 & FDR < 0.05 ~ "Upregulated",
        TRUE ~ "Non-significant"
      )
    )
  
  # Count regions by category for labels
  category_counts <- sig_regions %>%
    group_by(Category) %>%
    summarize(Count = n(), .groups = "drop")
  
  # Build legend labels with counts
  legend_labels <- category_counts %>%
    mutate(Label = paste(Category, "(", Count, ")")) %>%
    pull(Label)
  
  # Generate Volcano Plot
  volcano_plot <- ggplot(data = sig_regions, aes(x = logFC, y = -log10(FDR), color = Category)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Downregulated" = "darkgrey", 
                                  "Non-significant" = "lightblue", 
                                  "Upregulated" = "lightgrey"),
                       labels = legend_labels) +
    theme_classic() +
    ggtitle(paste(TF, phase, "Phase")) +
    theme(legend.position = "right") +
    xlim(c(-max(abs(results$logFC)), max(abs(results$logFC)))) +
    labs(x = "Log Fold Change", y = "-Log10 FDR", color = "Differential Expression")
  
  # Save the Volcano Plot
  ggsave(filename = paste0(file_path, TF, "_Volcano_Plot_", phase, ".pdf"), plot = volcano_plot)
  
  # Extract BED File
  results$chr <- sapply(strsplit(rownames(results), "[:-]"), `[`, 1)
  results$start <- as.numeric(sapply(strsplit(rownames(results), "[:-]"), `[`, 2))
  results$end <- as.numeric(sapply(strsplit(rownames(results), "[:-]"), `[`, 3))
  
  down_regions <- results %>%
    filter(logFC < -0.5 & FDR < 0.05) %>%
    select(chr, start, end) %>%
    mutate(category = "DOWN")
  
  up_regions <- results %>%
    filter(logFC > 0.5 & FDR < 0.05) %>%
    select(chr, start, end) %>%
    mutate(category = "UP")
  
  combined_regions <- bind_rows(down_regions, up_regions)
  bed_file <- paste0(file_path, TF, "_edgeR_DBRs_", phase, ".bed")
  write.table(combined_regions, file = bed_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # Region Count Plot
  lost <- nrow(down_regions)
  gained <- nrow(up_regions)
  count_data <- data.frame(
    Category = c("Lost Regions", "Gained Regions"),
    Count = c(lost, gained)
  )
  count_data <- count_data %>%
    mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))
  
  region_count_plot <- ggplot(count_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
    labs(title = paste("Number of", TF, "Regions in", phase, "Phase"),
         x = "", y = "Number of Regions", fill = "Differential Binding") +
    theme_minimal() +
    scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))
  
  # Save the Region Count Plot
  ggsave(filename = paste0(file_path, TF, "_Regions_", phase, ".pdf"), plot = region_count_plot)
}

# Analyze G1, S, and G2 Phases
analyze_phase(G1_merged, "G1", file_path, TF)
analyze_phase(S_merged, "S", file_path, TF)
analyze_phase(G2_merged, "G2", file_path, TF)

## Setup edgeR Parameters for summary table
analyze_phase_summary.table <- function(merged_data, phase, file_path, TF) {
  # Set up design and group
  group <- factor(c("ESC", "ESC", "EpiLC", "EpiLC"), levels = c("ESC", "EpiLC"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrast <- makeContrasts(reference_vs_other = "EpiLC - ESC", levels = design)
  
  # edgeR Differential Binding Analysis
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Fit the model
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract results
  results <- topTags(lrt, n = Inf)$table
  results$FDR <- p.adjust(results$PValue, method = "BH")
  return(results)
}

# Initialize an empty summary table
summary_table <- data.frame(
  Phase = character(),
  TF = character(),
  `ESC > EpiLC` = integer(),
  `ESC < EpiLC` = integer(),
  stringsAsFactors = FALSE
)

# Function to count regions and update the summary table
count_regions <- function(results, phase, tf, table) {
  down_count <- nrow(results %>% filter(logFC < -0.5 & FDR < 0.05))  # ESC > EpiLC
  up_count <- nrow(results %>% filter(logFC > 0.5 & FDR < 0.05))  # ESC < EpiLC
  
  # Add to summary table
  table <- rbind(
    table,
    data.frame(
      Phase = phase,
      TF = tf,
      `ESC > EpiLC` = down_count,
      `ESC < EpiLC` = up_count,
      stringsAsFactors = FALSE
    )
  )
  return(table)
}

# Analyze G1, S, and G2 Phases and generate the summary table
results_G1 <- analyze_phase_summary.table(G1_merged, "G1", file_path, TF)
summary_table <- count_regions(results_G1, "G1", TF, summary_table)

results_S <- analyze_phase_summary.table(S_merged, "S", file_path, TF)
summary_table <- count_regions(results_S, "S", TF, summary_table)

results_G2 <- analyze_phase_summary.table(G2_merged, "G2", file_path, TF)
summary_table <- count_regions(results_G2, "G2", TF, summary_table)

# Temporarily store the column names
column_names <- c("Phase", "TF", "ESC > EpiLC", "ESC < EpiLC")

# Write the column names manually and append the data
output_summary_file <- paste0(file_path, TF, "_Summary_Table.csv")
write(column_names, file = output_summary_file, ncolumns = length(column_names), sep = ",")
write.table(summary_table, file = output_summary_file, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE, append = TRUE)

# Print a success message
message("Summary table saved to: ", output_summary_file)
```
# Enrichment of ChIP-Atlas target genes
```R
# Read the TSV file into a data frame
df <- read.table("/Users/pascalhunold/Downloads/Nanog.1.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Extract the first column containing gene symbols
gene_symbols <- df$Target_genes

# Print or use gene_symbols as needed
print(gene_symbols)


# Path to your BED file
#mm39_bed <- "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/genomes/genes_mm39_both.sort.bed"
# Read BED file
#bed_data <- import.bed(mm39_bed)
# Convert to data frame for easier manipulation
#bed_df <- as.data.frame(bed_data)

# Filter BED data to keep only rows where gene symbol matches those in DESeq2 output
filtered_bed <- bed_df[bed_df$name %in% gene_symbols,]
# Select only the required columns in the order: chromosome, start, end, name, strand
output_bed <- filtered_bed[, c("seqnames", "start", "end", "name", "strand")]

# Write to a new BED file
write.table(output_bed, 
            file = "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/ChIP_Atlas/NANOG_ChIP_Atlas_mm10_Targets.bed", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


TF <- "NANOG"

# Read the data from the BED file into a data frame
peaks_df <- read.table(paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/peaks/", TF, "-G1_all_peaks.over59nt.sorted.bed"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
peaks_df$V1 <- paste0("chr", peaks_df$V1)
peaks <- GRanges(seqnames=peaks_df$V1, ranges=IRanges(start=peaks_df$V2, end=peaks_df$V3), strand=peaks_df$V4)

target_genes <- fread(paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/ChIP_Atlas/", TF, "_ChIP_Atlas_mm10_Targets.bed"), header = FALSE)
target_genes_gr <- GRanges(seqnames = target_genes$V1, 
                           ranges = IRanges(start = target_genes$V2, end = target_genes$V3))



# Read the known genes data
known_genes <- fread("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/genomes/genes_mm39_both.sort.bed", header = FALSE)

# Convert to GenomicRanges object
known_genes_gr <- GRanges(seqnames = known_genes$V1, 
                          ranges = IRanges(start = known_genes$V2, end = known_genes$V3))

# Extract gene symbols from known_genes
known_gene_symbols <- known_genes$V4

# Extract gene symbols from target_genes
target_gene_symbols <- target_genes$V4

# Find gene symbols in known_genes that are not in target_genes
non_target_gene_symbols <- setdiff(known_gene_symbols, target_gene_symbols)

# Randomly sample gene symbols from non_target_gene_symbols
# Determine the number of rows in target_genes
num_rows_target_genes <- nrow(target_genes)

# Randomly sample gene symbols from non_target_gene_symbols
random_gene_symbols <- sample(non_target_gene_symbols, num_rows_target_genes, replace = FALSE)
length(random_gene_symbols)
# Filter known_genes based on random_gene_symbols
random_genes <- known_genes[V4 %in% random_gene_symbols]

# Convert filtered data to GRanges
random_genes_gr <- GRanges(
  seqnames = random_genes$V1,
  ranges = IRanges(start = random_genes$V2, end = random_genes$V3),
  strand = random_genes$V6,
  names = random_genes$V4
)



overlap_target <- countOverlaps(peaks, target_genes_gr)
overlap_random <- countOverlaps(peaks, random_genes_gr)
overlap_target_random <- countOverlaps(target_genes_gr, random_genes_gr)
length(peaks)
length(target_genes_gr)
length(random_genes_gr)
cat("Overlap with target genes:", sum(overlap_target > 0), "\n")
cat("Overlap with randomly selected genes:", sum(overlap_random > 0), "\n")
cat("Overlap between target and randomly selected genes:", sum(overlap_target_random > 0), "\n")
#library(ggplot2)

# Create a data frame for plotting
plot_data <- data.frame(
  Category = c("ChIP Atlas Target Genes", "Random Genes"),
  Count = c(sum(overlap_target > 0), sum(overlap_random > 0))
)

# Plot
gg <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Overlap between Peaks and Genes -", TF),
       x = "",
       y = "Peaks in Genes [#]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Random Genes" = "#e3e3e3", "ChIP Atlas Target Genes" = "#666666")) +
  guides(fill = FALSE)
gg
# Save the plot as a PDF
ggsave(filename = paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/figures/", TF, "_ChIPAtlas_enrichment.pdf"), plot = gg)
```
# MEME-ChIP mm39
## Generate fasta Files
```bash
module load bedtools/2.29.2

for f in *.over59nt.sorted.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "module load bedtools/2.31.0 && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa -bed $f -fo ./${f%.conservedPeaks_MACS2_real.bed}.fasta"
done
```
## Generate fasta Files mm10
```bash
for f in *master_peaks.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa -bed $f -fo ./${f%.bed}.mm10.fasta"
done
```
## Generate fasta Files mm10 peaks_GSE224292_mESC_CUTnRUN
```bash
for f in *_peaks.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa -bed $f -fo ./${f%.bed}.mm10.fasta"
done
```
## Generate fasta Files mm10 peaks_GSE224292_mESC_CUTnRUN_no.control
```bash
for f in *peaks_no.control.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa -bed $f -fo ./${f%.bed}.mm10.fasta"
done
```
## Generate fasta Files mm10 GSM4291125_mESC_ChIPseq_peaks
```bash
for f in *_peaks.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa -bed $f -fo ./${f%.bed}.mm10.fasta"
done
```
## Generate fasta Files mm10 GSM4291125_mESC_ChIPseq_peaks_no_control
```bash
for f in *_peaks_no.control.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa -bed $f -fo ./${f%.bed}.mm10.fasta"
done
```
# TOBIAS Analysis
```bash
mkdir ATACorrect_ATAC_EpiLC

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# load own modules
module load use.own

# load TOBIAS
module load pypack/tobias

# run TOBIAS
tobias ATACorrect --bam EpiLC-ATAC.merged.sorted.bam --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir ./ATACorrect_ATAC_EpiLC --cores 8


# deactivate conda env
conda deactivate

---
mkdir ATACorrect_ATAC_ESC

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# load own modules
module load use.own

# load TOBIAS
module load pypack/tobias

# run TOBIAS
tobias ATACorrect --bam ESC-ATAC.merged.sorted.bam --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir ./ATACorrect_ATAC_ESC --cores 8


# deactivate conda env
conda deactivate


#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module purge
module load miniconda/py38_4.9.2
module load use.own
module load pypack/tobias
for f in *_corrected.bw
do
tobias FootprintScores --signal $f --regions /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --output ${f%%_corrected.bw}_footprints.bw --cores 8
done

mkdir BINDetect_output

#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module purge
module load miniconda/py38_4.9.2
module load use.own
module load pypack/tobias
tobias BINDetect --motifs /projects/ag-haensel/CUT_Tag/Gisela/tobias/motifs/JASPAR_all_motifs/all_motifs.jaspar --signals ESC-ATAC.merged.sorted_footprints.bw EpiLC-ATAC.merged.sorted_footprints.bw --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir BINDetect_output --cond_names ESC EpiLC --cores 8
```














