# Figure_3C-D
## snDynaTag_ESC_EpiLC_DynaTag_mm10_log2_sorted_and_coverage_in.snDynaTag_summits_filtered_by_merged_peaks_and_extended_by_half_avg_peak_size_filtered_1000.score_v1.sh
```bash
#!/bin/bash

# Directories
PEAKS_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/snDynaTag_GP_ESC_EpiLC_mm10_peaks/peaks/summits_filtered_1000_score_peaks_filtered_extended_by_avg_peak_size"
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/snDynaTag_GP_ESC_EpiLC_mm10_bigwig"
LOG2_BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/snDynaTag_GP_ESC_EpiLC_mm10_bigwig/log2_ratio"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/snDynaTag_ESC_EpiLC_DynaTag_mm10_log2_sorted_and_coverage_in.snDynaTag_summits_filtered_by_merged_peaks_and_extended_by_half_avg_peak_size_filtered_1000.score_v1"

# Mapping of file names to TF names
declare -A TF_MAP=(
    ["Myc"]="MYC"
    ["Nanog"]="NANOG"
    ["Pou5f1"]="OCT4"
    ["Yap1"]="YAP1"
)

# Mapping of zMin and zMax values for each TF
declare -A ZRANGE_MAP=(
    ["MYC"]="-2.6 2.6"
    ["NANOG"]="-1.8 1.8"
    ["OCT4"]="-1.8 1.8"
    ["YAP1"]="-2.1 2.1"
)

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Python script to filter matrix
FILTER_SCRIPT=$(cat <<'EOF'
import pandas as pd
import numpy as np
import sys

# Load matrix
matrix_file = sys.argv[1]
output_combined = sys.argv[2]

# Read matrix
data = pd.read_csv(matrix_file, sep="\t", skiprows=1, header=None)
columns = ['chr', 'start', 'end'] + [f'bin_{i}' for i in range(1, data.shape[1] - 3 + 1)]
data.columns = columns

# Replace 0s with NaNs to ignore in calculations
data.iloc[:, 3:] = data.iloc[:, 3:].replace(0, np.nan)

# Calculate mean log2 ratio (select numeric columns only)
data["mean_log2"] = data.iloc[:, 3:].select_dtypes(include=[np.number]).mean(axis=1, skipna=True)

# Filter regions
positive_regions = data[data["mean_log2"] > 0.5][["chr", "start", "end", "mean_log2"]]
negative_regions = data[data["mean_log2"] < -0.5][["chr", "start", "end", "mean_log2"]]

# Combine, drop duplicates, and save
combined_regions = pd.concat([positive_regions, negative_regions]).drop_duplicates().sort_values(by=["chr", "start"])
combined_regions.to_csv(output_combined, sep="\t", index=False, header=False)
EOF
)

# Write the Python filter script to the output directory
echo "$FILTER_SCRIPT" > "${OUTPUT_DIR}/filter_script.py"

# Array to store SLURM job IDs for filtering jobs
declare -a FILTER_JOB_IDS=()

# Initialize the region counts file
echo -e "TF\tPositive_Log2Regions\tNegative_Log2Regions" > "${OUTPUT_DIR}/region_counts.tsv"

# Function to submit job for each TF
submit_job() {
    local tf_key=$1
    local tf_name=${TF_MAP[$tf_key]}
    local peaks_file="${PEAKS_DIR}/${tf_name}_snDynaTag_extended_summits.bed"
    local matrix_file="${OUTPUT_DIR}/${tf_name}_log2_matrix.gz"
    local filtered_region_file="${OUTPUT_DIR}/${tf_name}_filtered_regions.bed"
    local filtered_matrix_file="${OUTPUT_DIR}/${tf_name}_filtered_matrix.gz"
    local sorted_regions_file="${OUTPUT_DIR}/${tf_name}_log2_sorted_regions.bed"
    local esc_bigwig="${BIGWIG_DIR}/merged_sample_${tf_name}_ESC_combined_mm10_cpm.bw"
    local epilc_bigwig="${BIGWIG_DIR}/merged_sample_${tf_name}_EpiLC_combined_mm10_cpm.bw"
    local log2_bigwig="${LOG2_BIGWIG_DIR}/log2_ratio_${tf_name}_ESC_vs_EpiLC_mm10.bw"

    local log2_heatmap_file="${OUTPUT_DIR}/${tf_name}_log2ratio_heatmap.pdf"
    local log2_log_file="${OUTPUT_DIR}/${tf_name}_log2_heatmap.log"
    local coverage_matrix_file="${OUTPUT_DIR}/${tf_name}_coverage_matrix.gz"
    local coverage_heatmap_file="${OUTPUT_DIR}/${tf_name}_coverage_heatmap.pdf"
    local coverage_log_file="${OUTPUT_DIR}/${tf_name}_coverage_heatmap.log"
    local plot_title="${tf_name} occupancy (snDynaTag) in ESC vs EpiLC ext.summits filtered (1000)"

    # Check if the input peaks file exists
    if [ ! -f "$peaks_file" ]; then
        echo "Error: Peaks file $peaks_file not found. Exiting."
        exit 1
    fi

    # Retrieve zMin and zMax values
    local z_range=${ZRANGE_MAP[$tf_name]}
    local z_min=$(echo $z_range | cut -d' ' -f1)
    local z_max=$(echo $z_range | cut -d' ' -f2)

    # Submit computeMatrix job for log2 ratio heatmap
    compute_matrix_job_id=$(sbatch --parsable --time=01:00:00 --mem=32gb --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    computeMatrix reference-point \
        --regionsFileName \"$peaks_file\" \
        --scoreFileName \"$log2_bigwig\" \
        --outFileName \"$matrix_file\" \
        --referencePoint center \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --binSize 50 --averageTypeBins mean --missingDataAsZero")

    # Run filtering after computeMatrix completes
    filter_job_id=$(sbatch --parsable --dependency=afterok:$compute_matrix_job_id --mem=64gb --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    python ${OUTPUT_DIR}/filter_script.py \"$matrix_file\" \"$filtered_region_file\" && \
    computeMatrix reference-point \
        --regionsFileName \"$filtered_region_file\" \
        --scoreFileName \"$log2_bigwig\" \
        --outFileName \"$filtered_matrix_file\" \
        --referencePoint center \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --averageTypeBins mean --missingDataAsZero --binSize 50 && \
    plotHeatmap \
        --matrixFile \"$filtered_matrix_file\" \
        --outFileName \"$log2_heatmap_file\" \
        --plotFileFormat pdf \
        --dpi 200 \
        --colorMap \"bwr\" \
        --samplesLabel \"\" \
        --legendLocation best \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow \"heatmap and colorbar\" \
        --sortRegions descend \
        --outFileSortedRegions \"$sorted_regions_file\" \
        --zMin \"$z_min\" \
        --zMax \"$z_max\" \
        --startLabel \"Upstream\" \
        --endLabel \"Downstream\" \
        --refPointLabel \"Summit center\" \
        --plotTitle \"$plot_title\" > \"$log2_log_file\" 2>&1 && \
    computeMatrix reference-point \
        --regionsFileName \"$sorted_regions_file\" \
        --scoreFileName \"$esc_bigwig\" \"$epilc_bigwig\" \
        --outFileName \"$coverage_matrix_file\" \
        --samplesLabel \"ESC\" \"EpiLC\" \
        --numberOfProcessors 8 \
        --referencePoint center \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --averageTypeBins \"mean\" \
        --missingDataAsZero \
        --binSize 50 && \
    plotHeatmap \
        --matrixFile \"$coverage_matrix_file\" \
        --outFileName \"$coverage_heatmap_file\" \
        --plotFileFormat pdf \
        --dpi 200 \
        --colorMap \"viridis\" \
        --samplesLabel \"ESC\" \"EpiLC\" \
        --legendLocation best \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow \"heatmap and colorbar\" \
        --sortRegions no \
        --startLabel \"Upstream\" \
        --endLabel \"Downstream\" \
        --refPointLabel \"Summit center\" \
        --plotTitle \"$plot_title\" > \"$coverage_log_file\" 2>&1 && \
    positive_count=\$(awk '\$4 > 0.5' \"$filtered_region_file\" | wc -l) && \
    negative_count=\$(awk '\$4 < -0.5' \"$filtered_region_file\" | wc -l) && \
    echo -e \"${tf_name}\t\${positive_count}\t\${negative_count}\" >> \"${OUTPUT_DIR}/region_counts.tsv\"")

    FILTER_JOB_IDS+=($filter_job_id)
}

# Submit jobs for all TFs
for tf_key in "${!TF_MAP[@]}"; do
    submit_job "$tf_key"
done

# Wait for all filtering jobs to complete before finishing
for job_id in "${FILTER_JOB_IDS[@]}"; do
    scontrol wait $job_id
done

echo "region_counts.tsv created successfully in ${OUTPUT_DIR}."
```
## filtered_1000_score_peaks_that_overlap_with_diff.log2ratio.ext.summits.bed.sh
```bash
#!/bin/bash

# Directories containing summit and peak files
SUMMIT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/snDynaTag_ESC_EpiLC_DynaTag_mm10_log2_sorted_and_coverage_in.snDynaTag_summits_filtered_by_merged_peaks_and_extended_by_half_avg_peak_size_filtered_1000.score_v1"
PEAK_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/snDynaTag_GP_ESC_EpiLC_mm10_peaks/peaks/snDynaTag_mm10_ESC_EpiLC_filtered_1000_score_peaks_per_TF_merged"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/snDynaTag_ESC_EpiLC_DynaTag_mm10_log2_sorted_and_coverage_in.snDynaTag_summits_filtered_by_merged_peaks_and_extended_by_half_avg_peak_size_filtered_1000.score_v1"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Transcription factors to process
TFS=("MYC" "NANOG" "OCT4")

# Activate the conda environment with bedtools
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Iterate over each transcription factor
for TF in "${TFS[@]}"; do
  SUMMIT_FILE="$SUMMIT_DIR/${TF}_filtered_regions.bed"
  PEAK_FILE="$PEAK_DIR/${TF}_snDynaTag_mm10_ESC_EpiLC_filtered_1000_score_merged_peaks.bed"
  OUTPUT_FILE="$OUTPUT_DIR/${TF}_filtered_1000_score_peaks_that_overlap_with_diff.log2ratio.ext.summits.bed"

  # Check if both the summit and peak files exist
  if [[ -f "$SUMMIT_FILE" && -f "$PEAK_FILE" ]]; then
    echo "Processing $TF..."

    # Use bedtools intersect to find overlapping peaks and only keep entries from -a
    bedtools intersect -a "$PEAK_FILE" -b "$SUMMIT_FILE" -u > "$OUTPUT_FILE"

    echo "Output written to $OUTPUT_FILE"
  else
    echo "Error: Summit or peak file for $TF not found. Skipping."

  fi
done

echo "Processing complete."
```
