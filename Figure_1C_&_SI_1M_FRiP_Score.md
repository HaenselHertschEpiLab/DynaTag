# FRiP Score
## Calculate FRiP Scores
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
## Plotting 
```R
# Read the file
data <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/FRiP_scores.txt", header = FALSE, stringsAsFactors = FALSE)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Extracting group names from the 'Sample' column
data$Group <- gsub("^(.*?)\\-.+?$", "\\1", data$V2)

# Plotting the box plot
ggplot(data, aes(x = Group, y = as.numeric(V6))) +
  geom_boxplot() +
  labs(x = "Group", y = "FRiP Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0, 1)
```
