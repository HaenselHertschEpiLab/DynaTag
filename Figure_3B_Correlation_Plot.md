# Correlation snDynaTag vs bulk ESC & EpiLC
```bash
#!/bin/bash

# Conda environment activation
CONDA_ENV="/projects/ag-haensel/tools/.conda/envs/abc-model-env"

# Paths to BAM files
SN_BAM_PATH="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/snDynaTag_GP_ESC_EpiLC_mm10_bam"
BULK_BAM_PATH="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"

# Define transcription factors and cell states
TFs=("NANOG" "MYC" "OCT4")
CELL_STATES=("ESC" "EpiLC")

# Output directories for results
CORRELATION_OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotCorrelation"
mkdir -p $CORRELATION_OUTPUT_DIR

# Bin size
BIN_SIZE=50

# Loop through each transcription factor and cell state
for TF in "${TFs[@]}"; do
  for STATE in "${CELL_STATES[@]}"; do

    # Define file paths for snDynaTag and bulkDynaTag BAM files
    SN_BAM="$SN_BAM_PATH/${TF}_${STATE}_merged_mm10.sorted.bam"

    if [[ $STATE == "EpiLC" ]]; then
      BULK_BAM="$BULK_BAM_PATH/EpiLC-d2-${TF}-G1.mm10.merged.sorted.bam"
    else
      BULK_BAM="$BULK_BAM_PATH/${STATE}-${TF}-G1.mm10.merged.sorted.bam"
    fi

    # Ensure bulkDynaTag BAM file exists and contains 'G1' in the filename
    if [[ -f $BULK_BAM && $BULK_BAM == *"G1"* ]]; then
      # Define output file names for bins jobs
      NPZ_OUTPUT_BINS="$CORRELATION_OUTPUT_DIR/${TF}_${STATE}_sn_vs_bulk_bins.npz"
      PDF_OUTPUT_BINS="$CORRELATION_OUTPUT_DIR/${TF}_${STATE}_sn_vs_bulk_bins_correlation.pdf"
      PDF_OUTPUT_BINS_NO_OUTLIERS="$CORRELATION_OUTPUT_DIR/${TF}_${STATE}_sn_vs_bulk_bins_correlation_no_outliers.pdf"

      # Define sample labels
      LABELS="${TF}_${STATE}_snDynaTag ${TF}_${STATE}_G1_bulkDynaTag"

      # Submit multiBamSummary job (bins) and capture job ID
      echo "Submitting bins job for ${TF} in ${STATE}..."
      JOB_ID_BINS=$(sbatch --parsable --mem=16G -J mBaSum_bins_${TF}_${STATE} --cpus-per-task=8 --wrap \
        "conda activate $CONDA_ENV && multiBamSummary bins -bs $BIN_SIZE -b $SN_BAM $BULK_BAM -p 8 -o $NPZ_OUTPUT_BINS")
      echo "JOB_ID_BINS: $JOB_ID_BINS"

      if [[ -z "$JOB_ID_BINS" ]]; then
        echo "Error: Failed to capture job ID for bins job. Skipping plotCorrelation jobs for ${TF} in ${STATE}."
        continue
      fi

      # Submit plotCorrelation job (with outliers) for bins
      echo "Submitting plotCorrelation (with outliers) job for bins ${TF} in ${STATE}..."
      sbatch --mem=64G --cpus-per-task=8 -J pCorr_bins_${TF}_${STATE} --dependency=afterok:${JOB_ID_BINS} --wrap \
        "conda activate $CONDA_ENV && \
        plotCorrelation -in $NPZ_OUTPUT_BINS --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p --labels $LABELS -o $PDF_OUTPUT_BINS"

      # Submit plotCorrelation job (without outliers) for bins
      echo "Submitting plotCorrelation (without outliers) job for bins ${TF} in ${STATE}..."
      sbatch --mem=64G --cpus-per-task=8 -J pCorr_bins_noOut_${TF}_${STATE} --dependency=afterok:${JOB_ID_BINS} --wrap \
        "conda activate $CONDA_ENV && \
        plotCorrelation -in $NPZ_OUTPUT_BINS --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p --removeOutliers --labels $LABELS -o $PDF_OUTPUT_BINS_NO_OUTLIERS"
    else
      echo "Bulk BAM file for ${TF} in ${STATE} state is missing or does not contain 'G1' in the filename. Skipping..."
    fi

  done
done
```
