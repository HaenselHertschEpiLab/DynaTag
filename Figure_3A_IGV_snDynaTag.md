## Merge replicates of aggregated BAM files together and bigwig generation mm10
```bash
#!/bin/bash

# Directory containing the BAM files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/snDynaTag_GP_ESC_EpiLC_mm10_bam"

# Output directory for merged BAM files
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/snDynaTag_GP_ESC_EpiLC_mm10_bam"
mkdir -p "$OUTPUT_DIR"

# Samples to merge
declare -A samples
samples=(
  ["MYC_ESC"]="merged_sample_MYC_ESC_1_mm10.sorted.bam merged_sample_MYC_ESC_2_mm10.sorted.bam"
  ["MYC_EpiLC"]="merged_sample_MYC_EpiLC_1_mm10.sorted.bam merged_sample_MYC_EpiLC_2_mm10.sorted.bam"
  ["NANOG_ESC"]="merged_sample_NANOG_ESC_1_mm10.sorted.bam merged_sample_NANOG_ESC_2_mm10.sorted.bam"
  ["NANOG_EpiLC"]="merged_sample_NANOG_EpiLC_1_mm10.sorted.bam merged_sample_NANOG_EpiLC_2_mm10.sorted.bam"
  ["OCT4_ESC"]="merged_sample_OCT4_ESC_1_mm10.sorted.bam merged_sample_OCT4_ESC_2_mm10.sorted.bam"
  ["OCT4_EpiLC"]="merged_sample_OCT4_EpiLC_1_mm10.sorted.bam merged_sample_OCT4_EpiLC_2_mm10.sorted.bam"
  ["YAP1_ESC"]="merged_sample_YAP1_ESC_1_mm10.sorted.bam merged_sample_YAP1_ESC_2_mm10.sorted.bam"
  ["YAP1_EpiLC"]="merged_sample_YAP1_EpiLC_1_mm10.sorted.bam merged_sample_YAP1_EpiLC_2_mm10.sorted.bam"
)

# Activate software
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Loop through samples and merge replicates
for sample in "${!samples[@]}"; do
  input_files=""
  for file in ${samples[$sample]}; do
    input_files+="$BAM_DIR/$file "
  done

  # Output file
  merged_file="$OUTPUT_DIR/${sample}_merged_mm10.bam"
  sorted_file="$OUTPUT_DIR/${sample}_merged_mm10.sorted.bam"

  echo "Merging files for $sample into $merged_file"
  samtools merge -@ 4 "$merged_file" $input_files && \
  echo "Sorting $merged_file into $sorted_file" && \
  samtools sort -@ 4 -o "$sorted_file" "$merged_file" && \
  echo "Indexing $sorted_file" && \
  samtools index "$sorted_file" && \
  echo "Removing intermediate file $merged_file" && \
  rm "$merged_file"

done

# Notify that merging is complete
echo "All BAM files have been merged, sorted, and indexed. Starting BigWig file generation..."

# Generate BigWig files
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/snDynaTag_GP_ESC_EpiLC_mm10_bigwig"
mkdir -p "$BIGWIG_DIR"

for f1 in "$OUTPUT_DIR"/*_merged_mm10.sorted.bam; do
  echo "Generating BigWig file for $f1"
  bamCoverage -p 8 -b "$f1" -bs 10 --skipNAs --centerReads --normalizeUsing CPM -of bigwig -o "$BIGWIG_DIR/$(basename ${f1%%.sorted.bam})_cpm.bw"
done

# Count reads and save to text file
READ_COUNT_FILE="$OUTPUT_DIR/read_counts.txt"
echo -e "File\tReadCount" > "$READ_COUNT_FILE"
for f1 in "$OUTPUT_DIR"/*_merged_mm10.sorted.bam; do
  read_count=$(samtools view -c "$f1")
  echo -e "$(basename "$f1")\t$read_count" >> "$READ_COUNT_FILE"
  echo "Counted $read_count reads in $f1"
done

echo "All samples merged, sorted, indexed, BigWig files generated, and read counts saved to $READ_COUNT_FILE successfully."
```
