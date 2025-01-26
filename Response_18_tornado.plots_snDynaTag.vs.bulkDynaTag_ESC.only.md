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








