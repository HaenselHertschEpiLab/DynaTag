# Figure_2A-C_DOA_ESC_vs_EpiLC_DynaTag.peaks_mm10
```bash
cat Generate_Count_Matrices_mm10.sh
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
```

```bash
cat Generate_Count_Matrices_mm10_not.norm.sh
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

```bash
cat rename_counts.txt.sh
for file in *_mm10_counts.txt; do
  # Extract the base name before "_mm10"
  new_name=$(echo "$file" | sed -E 's/_S[0-9]+_[0-9]+_mm10_counts.txt/_counts.txt/')
  # Rename the file
  mv "$file" "$new_name"
done
```

```bash
# R analysis edgeR DynaTag peaks example for NANOG
# RHH reproduces Pascal's Differential Binding Analysis mm10
## Load Packages
library(edgeR)
library(ggplot2)
library(dplyr)

## Load Data and Assemble Count Matrices
rm(list = ls())

# Define the file path
file_path <- "/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices_mm10/"

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

sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4   ggplot2_3.5.1 edgeR_4.2.2   limma_3.60.6 

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5      cli_3.6.3        rlang_1.1.4      generics_0.1.3   glue_1.8.0      
 [6] labeling_0.4.3   statmod_1.5.0    colorspace_2.1-1 locfit_1.5-9.10  scales_1.3.0    
[11] grid_4.4.1       munsell_0.5.1    tibble_3.2.1     lifecycle_1.0.4  compiler_4.4.1  
[16] Rcpp_1.0.14      pkgconfig_2.0.3  farver_2.1.2     lattice_0.22-6   R6_2.5.1        
[21] tidyselect_1.2.1 pillar_1.10.1    splines_4.4.1    magrittr_2.0.3   tools_4.4.1     
[26] withr_3.0.2      gtable_0.3.6
```

```bash
cat bulkDynaTag_ESC_EpiLC_DynaTag_mm10_log2_analysis_keep.sort_then_binding_in_edgeR_mm10_summits.of.peaks_v2.sh
#!/bin/bash

# Directories
PEAKS_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/edgeR/mm10_edgeR"
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_bigwig"
LOG2_BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_bigwig/log2_ratio_bulk"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_log2_analysis_keep.sort_then_binding_in_edgeR_mm10_summits.of.peaks_v2"

# Mapping of file names to TF names
declare -A TF_MAP=(
    ["MYC"]="MYC"
    ["NANOG"]="NANOG"
    ["OCT4"]="OCT4"
    ["SOX2"]="SOX2"
    ["YAP1"]="YAP1"
)

# Mapping of zMin and zMax values for each TF
declare -A ZRANGE_MAP=(
    ["MYC"]="-2.2 2.2"
    ["NANOG"]="-1.5 1.5"
    ["OCT4"]="-1.5 1.5"
    ["SOX2"]="-1.2 1.2"
    ["YAP1"]="-2 2"
)

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Phases to analyze
PHASES=("G1" "S" "G2")

# Function to submit job for each TF and phase
submit_job() {
    local tf_name=$1
    local phase=$2

    # Construct file paths
    local region_file="${PEAKS_DIR}/${tf_name}_edgeR_DBRs_${phase}.summits.bed"
    local matrix_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_matrix.gz"
    local sorted_regions_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_sorted_regions.bed"
    local esc_bigwig="${BIGWIG_DIR}/ESC-${tf_name}-${phase}.mm10.merged_cpm.bw"
    local epilc_bigwig="${BIGWIG_DIR}/EpiLC-d2-${tf_name}-${phase}.mm10.merged_cpm.bw"

    local log2_heatmap_file="${OUTPUT_DIR}/${tf_name}_log2ratio_heatmap_${phase}.pdf"
    local log2_log_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_heatmap.log"
    local coverage_matrix_file="${OUTPUT_DIR}/${phase}_${tf_name}_coverage_matrix.gz"
    local coverage_heatmap_file="${OUTPUT_DIR}/${tf_name}_coverage_heatmap_${phase}.pdf"
    local coverage_log_file="${OUTPUT_DIR}/${phase}_${tf_name}_coverage_heatmap.log"
    local plot_title="${tf_name} occupancy (bulkDynaTag) in ESC EpiLC ${phase} diff. edgeR summits"

    # Retrieve zMin and zMax values
    local z_range=${ZRANGE_MAP[$tf_name]}
    local z_min=$(echo $z_range | cut -d' ' -f1)
    local z_max=$(echo $z_range | cut -d' ' -f2)

    # Submit computeMatrix and plotHeatmap jobs
    sbatch --time=01:00:00 --mem=32gb --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    computeMatrix reference-point \
        --regionsFileName \"$region_file\" \
        --scoreFileName \"$LOG2_BIGWIG_DIR/log2_ratio_ESC-${tf_name}-${phase}.mm10.merged_cpm_vs_EpiLC-d2-${tf_name}-${phase}.mm10.merged_cpm.bw\" \
        --outFileName \"$matrix_file\" \
        --referencePoint center \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --binSize 50 --averageTypeBins mean --missingDataAsZero && \
    plotHeatmap \
        --matrixFile \"$matrix_file\" \
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
        --refPointLabel \"Peak center\" \
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
        --refPointLabel \"Peak center\" \
        --plotTitle \"$plot_title\" > \"$coverage_log_file\" 2>&1"
}

# Submit jobs for all TFs and phases
for phase in "${PHASES[@]}"; do
    for tf_name in "${!TF_MAP[@]}"; do
        submit_job "$tf_name" "$phase"
    done
done
```
# Figure_2A-C_DOA_ESC_vs_EpiLC_otherTFs_DynaTag.peaks_mm10
```bash
Same procedure for other TFs MYC, YAP1, OCT4, SOX2
```
# Figure_SI_3B-D_DOA_ESC_vs_EpiLC_MYC_ChIPseq.peaks_mm10
```bash
## Generation of differential occupied regions by edgeR mm10 --> PEAKS_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/edgeR/mm10_edgeR_ChIPseq_peaks"
cat Generate_Count_Matrices_ChIP.Atlas.peaks_mm10.sh
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ChIP.Atlas_TFs/3col.peaks_ChIP.Atlas_TFs"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_ChIP.Atlas.peaks_mm10"

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/merged_Oth.PSC.05.${epitope}.AllCell.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*_mm10_norm_clean.sort.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph"
                    awk -v OFS='\t' '{print $1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%_norm_clean.sort.bam})_counts.txt"
                done
            fi
        done
    done
done
```

```bash
cat Generate_Count_Matrices_ChIP.Atlas.peaks_mm10_not.norm.sh
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_ChIP.Atlas_TFs/3col.peaks_ChIP.Atlas_TFs"
bam_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"
bedgraph_path="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bedgraph/bedgraph_ChIP.Atlas.peaks_mm10"
for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/merged_Oth.PSC.05.${epitope}.AllCell.bed"
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

```bash
cat rename_counts.txt.sh
for file in *_mm10_counts.txt; do
  # Extract the base name before "_mm10"
  new_name=$(echo "$file" | sed -E 's/_S[0-9]+_[0-9]+_mm10_counts.txt/_counts.txt/')
  # Rename the file
  mv "$file" "$new_name"
done
```

```bash
# R analysis edgeR known ChIP-seq peaks
## Load Packages
library(edgeR)
library(ggplot2)
library(dplyr)

## Load Data and Assemble Count Matrices
rm(list = ls())

# Define the file path
file_path <- "/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices_ChIP.Atlas_mm10/"

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

R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4   ggplot2_3.5.1 edgeR_4.2.2   limma_3.60.6 

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5      cli_3.6.3        rlang_1.1.4      generics_0.1.3   glue_1.8.0       labeling_0.4.3   statmod_1.5.0   
 [8] colorspace_2.1-1 locfit_1.5-9.10  scales_1.3.0     grid_4.4.1       munsell_0.5.1    tibble_3.2.1     lifecycle_1.0.4 
[15] compiler_4.4.1   Rcpp_1.0.14      pkgconfig_2.0.3  farver_2.1.2     lattice_0.22-6   R6_2.5.1         tidyselect_1.2.1
[22] pillar_1.10.1    splines_4.4.1    magrittr_2.0.3   tools_4.4.1      withr_3.0.2      gtable_0.3.6
```

```bash
cat bulkDynaTag_ESC_EpiLC_DynaTag_mm10_log2_analysis_keep.sort_then_binding_in_edgeR_mm10_ChIPseq.peaks_v2.sh
#!/bin/bash

# Directories
PEAKS_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/edgeR/mm10_edgeR_ChIPseq_peaks"
BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_bigwig"
LOG2_BIGWIG_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/bigwig/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_bigwig/log2_ratio_bulk"
OUTPUT_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/plotHeatmap/bulkDynaTag_ESC_EpiLC_DynaTag_mm10_log2_analysis_keep.sort_then_binding_in_edgeR_mm10_ChIPseq.peaks_v2"

# Mapping of file names to TF names
declare -A TF_MAP=(
    ["MYC"]="MYC"
    ["NANOG"]="NANOG"
    ["OCT4"]="OCT4"
    ["SOX2"]="SOX2"
    ["YAP1"]="YAP1"
)

# Mapping of zMin and zMax values for each TF
declare -A ZRANGE_MAP=(
    ["MYC"]="-3 3"
    ["NANOG"]="-1.5 1.5"
    ["OCT4"]="-1.6 1.6"
    ["SOX2"]="-1.1 1.1"
    ["YAP1"]="-2.1 2.1"
)

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Phases to analyze
PHASES=("G1" "S" "G2")

# Function to submit job for each TF and phase
submit_job() {
    local tf_name=$1
    local phase=$2

    # Construct file paths
    local region_file="${PEAKS_DIR}/${tf_name}_edgeR_DBRs_${phase}.bed"
    local matrix_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_matrix.gz"
    local sorted_regions_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_sorted_regions.bed"
    local esc_bigwig="${BIGWIG_DIR}/ESC-${tf_name}-${phase}.mm10.merged_cpm.bw"
    local epilc_bigwig="${BIGWIG_DIR}/EpiLC-d2-${tf_name}-${phase}.mm10.merged_cpm.bw"

    local log2_heatmap_file="${OUTPUT_DIR}/${tf_name}_log2ratio_heatmap_${phase}.pdf"
    local log2_log_file="${OUTPUT_DIR}/${phase}_${tf_name}_log2_heatmap.log"
    local coverage_matrix_file="${OUTPUT_DIR}/${phase}_${tf_name}_coverage_matrix.gz"
    local coverage_heatmap_file="${OUTPUT_DIR}/${tf_name}_coverage_heatmap_${phase}.pdf"
    local coverage_log_file="${OUTPUT_DIR}/${phase}_${tf_name}_coverage_heatmap.log"
    local plot_title="${tf_name} occupancy (bulkDynaTag) in ESC EpiLC ${phase} diff. edgeR ChIPseq peaks"

    # Retrieve zMin and zMax values
    local z_range=${ZRANGE_MAP[$tf_name]}
    local z_min=$(echo $z_range | cut -d' ' -f1)
    local z_max=$(echo $z_range | cut -d' ' -f2)

    # Submit computeMatrix and plotHeatmap jobs
    sbatch --time=01:00:00 --mem=32gb --cpus-per-task=8 --wrap "
    conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env && \
    computeMatrix reference-point \
        --regionsFileName \"$region_file\" \
        --scoreFileName \"$LOG2_BIGWIG_DIR/log2_ratio_ESC-${tf_name}-${phase}.mm10.merged_cpm_vs_EpiLC-d2-${tf_name}-${phase}.mm10.merged_cpm.bw\" \
        --outFileName \"$matrix_file\" \
        --referencePoint center \
        --beforeRegionStartLength 2500 \
        --afterRegionStartLength 2500 \
        --binSize 50 --averageTypeBins mean --missingDataAsZero && \
    plotHeatmap \
        --matrixFile \"$matrix_file\" \
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
        --startLabel \"Upstream\" \
        --endLabel \"Downstream\" \
        --zMin \"$z_min\" \
        --zMax \"$z_max\" \
        --refPointLabel \"Peak center\" \
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
        --refPointLabel \"Peak center\" \
        --plotTitle \"$plot_title\" > \"$coverage_log_file\" 2>&1"
}

# Submit jobs for all TFs and phases
for phase in "${PHASES[@]}"; do
    for tf_name in "${!TF_MAP[@]}"; do
        submit_job "$tf_name" "$phase"
    done
done
```
