# Response_18
## CHIP-ATLAS target gene processing mm10
```bash
mkdir peaks_consensus_target.genes_CHIP-ATLAS
nano process_multiple_targets.sh
#CHIP-ATLAS target genes obtained from https://chip-atlas.org/peak_browser

#!/bin/bash

# Define input files
coordinate_file="Genes_mm10_both.sort.bed"
target_gene_files=(
    "MYC_mESC_mm10_target.genes_10kb.txt"
    "Nanog_mESC_mm10_target.genes_10kb.txt"
    "OCT4_mESC_mm10_target.genes_10kb.txt"
    "SOX2_mESC_mm10_target.genes_10kb.txt"
    "YAP1_mESC_mm10_target.genes_10kb.txt"
)

# Output directory for results
output_dir="filtered_genes"
mkdir -p "$output_dir"

# Process each target gene file
for target_genes_file in "${target_gene_files[@]}"; do
    # Extract base name of the target gene file for naming outputs
    base_name=$(basename "$target_genes_file" .txt)

    # Extract target gene names in order (first column after the header)
    tail -n +2 "$target_genes_file" | cut -f1 > sorted_target_genes.txt

    # Use awk to filter the BED file based on target genes
    awk 'NR==FNR {genes[$1]; next} $4 in genes' sorted_target_genes.txt "$coordinate_file" > "$output_dir/${base_name}_filtered_genes.bed"

    # Optional: Sort filtered results to match the order of the target genes
    awk 'NR==FNR {order[$1]=NR; next} {if ($4 in order) print order[$4], $0}' sorted_target_genes.txt "$output_dir/${base_name}_filtered_genes.bed" | \
    sort -n | cut -d' ' -f2- > "$output_dir/${base_name}_ordered_filtered_genes.bed"

    echo "Processed $target_genes_file: Results saved in '$output_dir/${base_name}_ordered_filtered_genes.bed'"
done

# Clean up intermediate file
rm sorted_target_genes.txt
```
## Compute Matrix && plotHeatmap analysis of snDynaTag ESC signal at TSS of target genes. snDynaTag ESC experiments performed by PH.
```bash
nano computematrix_plotHeatmap_snDynaTag_ESC.only.sh

#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Directories
REGION_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_consensus_target.genes_CHIP-ATLAS/filtered_genes"
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/snDynaTag_PH_ESC.only_mm10_bigwig"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/snDynaTag_pseudobulk_plotHeatmaps"

# TF files
TFS=(
    "MYC_mESC_mm10_target.genes_10kb_ordered_filtered_genes.bed"
    "Nanog_mESC_mm10_target.genes_10kb_ordered_filtered_genes.bed"
    "OCT4_mESC_mm10_target.genes_10kb_ordered_filtered_genes.bed"
    "YAP1_mESC_mm10_target.genes_10kb_ordered_filtered_genes.bed"
)

# Function to process each TF
do_analysis() {
    local tf_name=$1
    local region_file=$2
    local bigwig_file=$3

    # Basename for outputs
    local base_name=$(basename "$region_file" .bed)

    # Output files
    local matrix_file_500="${OUTPUT_DIR}/${base_name}_matrix_500.gz"
    local matrix_file_2500="${OUTPUT_DIR}/${base_name}_matrix_2500.gz"
    local heatmap_file_500="${OUTPUT_DIR}/${base_name}_heatmap_500.pdf"
    local heatmap_file_2500="${OUTPUT_DIR}/${base_name}_heatmap_2500.pdf"

    # ComputeMatrix: +/- 500 bp
    echo "Running computeMatrix for +/- 500 bp: $tf_name"
    computeMatrix reference-point \
        --regionsFileName "$region_file" \
        --scoreFileName "$bigwig_file" \
        --outFileName "$matrix_file_500" \
        --samplesLabel "snDynaTag (Pseudobulk) $tf_name" \
        --numberOfProcessors 8 \
        --referencePoint TSS \
        --beforeRegionStartLength 500 \
        --afterRegionStartLength 500 \
        --sortRegions "keep" \
        --averageTypeBins "mean" \
        --missingDataAsZero \
        --binSize 50

    # PlotHeatmap: +/- 500 bp
    echo "Running plotHeatmap for +/- 500 bp: $tf_name"
    plotHeatmap --matrixFile "$matrix_file_500" \
        --outFileName "$heatmap_file_500" \
        --plotFileFormat "pdf" \
        --dpi "200" \
        --sortRegions "no" \
        --colorMap binary \
        --alpha "1.0" \
        --zMin 0 \
        --xAxisLabel "Distance" \
        --yAxisLabel "Target_Genes" \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow "heatmap and colorbar" \
        --startLabel "Upstream" \
        --endLabel "Downstream" \
        --refPointLabel "TSS" \
        --samplesLabel "snDynaTag (Pseudobulk) $tf_name" \
                --legendLocation "best" \
        --labelRotation "0"

    # ComputeMatrix: +/- 2500 bp
    echo "Running computeMatrix for +/- 2500 bp: $tf_name"
    computeMatrix reference-point \
        --regionsFileName "$region_file" \
        --scoreFileName "$bigwig_file" \
        --outFileName "$matrix_file_2500" \
        --samplesLabel "snDynaTag (Pseudobulk) $tf_name" \
        --numberOfProcessors 8 \
        --referencePoint TSS \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --sortRegions "keep" \
        --averageTypeBins "mean" \
        --missingDataAsZero \
        --binSize 50

    # PlotHeatmap: +/- 2500 bp
    echo "Running plotHeatmap for +/- 2500 bp: $tf_name"
    plotHeatmap --matrixFile "$matrix_file_2500" \
        --outFileName "$heatmap_file_2500" \
        --plotFileFormat "pdf" \
        --dpi "200" \
        --sortRegions "no" \
        --colorMap binary \
        --alpha "1.0" \
        --zMin 0 \
        --xAxisLabel "Distance" \
        --yAxisLabel "Target_Genes" \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow "heatmap and colorbar" \
        --startLabel "Upstream" \
        --endLabel "Downstream" \
        --refPointLabel "TSS" \
        --samplesLabel "snDynaTag (Pseudobulk) $tf_name" \
                --legendLocation "best" \
        --labelRotation "0"
}

# Main script
for tf_file in "${TFS[@]}"; do
    tf_name=$(basename "$tf_file" | cut -d'_' -f1 | tr '[:upper:]' '[:lower:]')
    bigwig_name="merged_sample_${tf_name^^}_mm10_cpm.bw"
    region_file="${REGION_DIR}/${tf_file}"
    bigwig_file="${BIGWIG_DIR}/${bigwig_name}"

    if [[ -f "$region_file" && -f "$bigwig_file" ]]; then
        do_analysis "$tf_name" "$region_file" "$bigwig_file"
    else
        echo "Missing files for $tf_name: $region_file or $bigwig_file"
    fi
done
```
## Compute Matrix && plotHeatmap analysis of snDynaTag ESC and EpiLC signal at TSS of target genes. snDynaTag experiments performed by GP and OvR.
```bash
nano computematrix_plotHeatmap_snDynaTag_GP_ESC_EpiLC_mm10.sh

#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rhaensel@uni-koeln.de

# Activate conda environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Directories
REGION_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_consensus_target.genes_CHIP-ATLAS/filtered_genes"
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/snDynaTag_GP_ESC_EpiLC_mm10_bigwig"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/snDynaTag_pseudobulk_GP_ESC_EpiLC_mm10_plotHeatmaps"

# TF files
TFS=(
    "MYC_ESC_1"
    "MYC_ESC_2"
    "MYC_EpiLC_1"
    "MYC_EpiLC_2"
    "NANOG_ESC_1"
    "NANOG_ESC_2"
    "NANOG_EpiLC_1"
    "NANOG_EpiLC_2"
    "OCT4_ESC_1"
    "OCT4_ESC_2"
    "OCT4_EpiLC_1"
    "OCT4_EpiLC_2"
    "YAP1_ESC_1"
    "YAP1_ESC_2"
    "YAP1_EpiLC_1"
    "YAP1_EpiLC_2"
)

# Function to process each TF and condition
do_analysis() {
    local tf_name=$1
    local region_file=$2
    local bigwig_file=$3

    # Basename for outputs, including TF name
    local base_name=$(basename "$region_file" .bed)_${tf_name}

    # Output files
    local matrix_file_500="${OUTPUT_DIR}/${base_name}_matrix_500.gz"
    local matrix_file_2500="${OUTPUT_DIR}/${base_name}_matrix_2500.gz"
    local heatmap_file_500="${OUTPUT_DIR}/${base_name}_heatmap_500.pdf"
    local heatmap_file_2500="${OUTPUT_DIR}/${base_name}_heatmap_2500.pdf"

    # Debugging logs
    echo "Processing TF: $tf_name"
    echo "Region file: $region_file"
    echo "BigWig file: $bigwig_file"
    echo "Output matrix (500bp): $matrix_file_500"
    echo "Output heatmap (500bp): $heatmap_file_500"

    # ComputeMatrix: +/- 500 bp
    echo "Running computeMatrix for +/- 500 bp..."
    computeMatrix reference-point \
        --regionsFileName "$region_file" \
        --scoreFileName "$bigwig_file" \
        --outFileName "$matrix_file_500" \
        --samplesLabel "snDynaTag (Pseudobulk) ${tf_name}" \
        --numberOfProcessors 8 \
        --referencePoint TSS \
        --beforeRegionStartLength 500 \
        --afterRegionStartLength 500 \
        --sortRegions "keep" \
        --averageTypeBins "mean" \
        --missingDataAsZero \
        --binSize 50 > "${OUTPUT_DIR}/computeMatrix_${tf_name}_500.log" 2>&1 || echo "Error running computeMatrix for ${tf_name} (500bp)"

    # PlotHeatmap: +/- 500 bp
    echo "Running plotHeatmap for +/- 500 bp..."
    plotHeatmap --matrixFile "$matrix_file_500" \
        --outFileName "$heatmap_file_500" \
        --plotFileFormat "pdf" \
        --dpi "200" \
        --sortRegions "no" \
        --colorMap binary \
        --alpha "1.0" \
        --zMin 0 \
        --zMax 2 \
        --xAxisLabel "Distance" \
        --yAxisLabel "Target_Genes" \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow "heatmap and colorbar" \
        --startLabel "Upstream" \
        --endLabel "Downstream" \
        --refPointLabel "TSS" \
        --samplesLabel "snDynaTag (Pseudobulk) ${tf_name}" \
        --legendLocation "best" \
        --labelRotation "0" > "${OUTPUT_DIR}/plotHeatmap_${tf_name}_500.log" 2>&1 || echo "Error running plotHeatmap for ${tf_name} (500bp)"

    # ComputeMatrix: +/- 2500 bp
    echo "Running computeMatrix for +/- 2500 bp..."
    computeMatrix reference-point \
        --regionsFileName "$region_file" \
        --scoreFileName "$bigwig_file" \
        --outFileName "$matrix_file_2500" \
        --samplesLabel "snDynaTag (Pseudobulk) ${tf_name}" \
        --numberOfProcessors 8 \
        --referencePoint TSS \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --sortRegions "keep" \
        --averageTypeBins "mean" \
        --missingDataAsZero \
        --binSize 50 > "${OUTPUT_DIR}/computeMatrix_${tf_name}_2500.log" 2>&1 || echo "Error running computeMatrix for ${tf_name} (2500bp)"

    # PlotHeatmap: +/- 2500 bp
    plotHeatmap --matrixFile "$matrix_file_2500" \
        --outFileName "$heatmap_file_2500" \
        --plotFileFormat "pdf" \
        --dpi "200" \
        --sortRegions "no" \
        --colorMap binary \
        --alpha "1.0" \
        --zMin 0 \
        --zMax 0.25 \
        --xAxisLabel "Distance" \
        --yAxisLabel "Target_Genes" \
        --heatmapWidth 7.5 \
        --heatmapHeight 7.5 \
        --whatToShow "heatmap and colorbar" \
        --startLabel "Upstream" \
        --endLabel "Downstream" \
        --refPointLabel "TSS" \
        --samplesLabel "snDynaTag (Pseudobulk) ${tf_name}" \
        --legendLocation "best" \
        --labelRotation "0" > "${OUTPUT_DIR}/plotHeatmap_${tf_name}_2500.log" 2>&1 || echo "Error running plotHeatmap for ${tf_name} (2500bp)"
}

# Main script
for tf_name in "${TFS[@]}"; do
    tf_base=$(echo "$tf_name" | cut -d'_' -f1) # Extract base name (e.g., YAP1)
    region_file="${REGION_DIR}/${tf_base}_mESC_mm10_target.genes_10kb_ordered_filtered_genes.bed"
    bigwig_file="${BIGWIG_DIR}/merged_sample_${tf_name}_mm10_cpm.bw"

    echo "Processing TF: $tf_name"
    echo "Looking for BigWig file: $bigwig_file"

    if [[ -f "$bigwig_file" ]]; then
        do_analysis "$tf_name" "$region_file" "$bigwig_file"
    else
        echo "No matching BigWig file found for $tf_name in $BIGWIG_DIR"
    fi
done
```
