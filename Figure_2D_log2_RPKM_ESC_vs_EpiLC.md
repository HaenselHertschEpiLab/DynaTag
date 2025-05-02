# Figure_2D_log2_RPKM_ESC_vs_EpiLC
```bash
# sbatch run_multicov.slurm.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

# Directory with BAM files
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"

# Directory with BED files
BED_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/edgeR/mm10_edgeR"

# List of BAM files (adjust the order if needed)
BAMS=(
  "$BAM_DIR/ESC-MYC-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-YAP1-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-NANOG-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-OCT4-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-SOX2-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-MYC-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-NANOG-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-OCT4-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-SOX2-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-YAP1-G1.mm10.merged.sorted.bam"
)

# Process only *_Gained.bed and *_Lost.bed files in BED_DIR.
for BED in "$BED_DIR"/*_Gained.bed "$BED_DIR"/*_Lost.bed; do
  if [ -f "$BED" ]; then
    OUTPUT="${BED%.bed}_multicov.txt"

    # Build the header: first three columns from the BED file, then one per BAM file.
    header="chr\tstart\tend"
    for bam in "${BAMS[@]}"; do
      # Get the base filename and remove everything from '.mm10' onward.
      fname=$(basename "$bam")
      shortname=${fname%%.mm10*}
      header="${header}\t${shortname}"
    done

    # Write the header and then append the bedtools multicov results.
    echo -e "$header" > "$OUTPUT"
    bedtools multicov -bams "${BAMS[@]}" -bed "$BED" >> "$OUTPUT"

    echo "Coverage for $BED saved to $OUTPUT"
  fi
done
```
```bash
# sbatch normalize_multicov.slurm.sh

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb

# Activate the required environment
conda activate /projects/ag-haensel/tools/.conda/envs/abc-model-env

########################################
# Define directories and BAM file list #
########################################

# Directory with BAM files (same as used before)
BAM_DIR="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/alignment/bam/bulk_DynaTag_CUTnTag_ESC_EpiLC_bam"

# List of BAM files (the order here should match the order used in your multicov files)
BAMS=(
  "$BAM_DIR/ESC-MYC-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-YAP1-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-NANOG-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-OCT4-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/ESC-SOX2-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-MYC-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-NANOG-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-OCT4-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-SOX2-G1.mm10.merged.sorted.bam"
  "$BAM_DIR/EpiLC-d2-YAP1-G1.mm10.merged.sorted.bam"
)

########################################
# Compute library sizes for each BAM   #
########################################

# We build a comma-separated string of library sizes (total read count per BAM)
libSizes=""
for bam in "${BAMS[@]}"; do
  # Get a short name by stripping off the suffix starting at ".mm10"
  fname=$(basename "$bam")
  shortname=${fname%%.mm10*}
  # Use samtools to count the number of reads (adjust the options if needed)
  count=$(samtools view -c "$bam")
  if [ -z "$libSizes" ]; then
    libSizes="$count"
  else
    libSizes="${libSizes},${count}"
  fi
done

echo "Library sizes for BAM files (in order): $libSizes"

########################################
# Process each multicov file to       #
# compute normalized values            #
########################################

# Loop over all *_multicov.txt files produced by the previous script.
# These files are assumed to have a header like:
#   chr    start    end    ESC-MYC-G1    ESC-YAP1    ... (short names)
for file in *_multicov.txt; do
  if [ -f "$file" ]; then
    out_file="${file%.txt}_normalized.txt"
    echo "Processing $file -> $out_file"
    
    # Process with awk.
    # The awk script:
    #  - Splits the passed-in library sizes into an array (lib[1] corresponds to column 4, etc.)
    #  - For each row (after the header), calculates the peak length (end - start)
    #  - Then, for each BAM column, calculates:
    #       normalized = (raw_count * 1e9) / (peak_length * library_size)
    #  - The header line is re-written to prepend "norm_" to each BAM column name.
    awk -v libSizes="$libSizes" 'BEGIN {
      OFS="\t";
      # Split the comma-separated library sizes into array "lib"
      n = split(libSizes, lib, ",");
    }
    NR==1 {
      # Process header: keep first three columns as is, and for columns 4...NF add a prefix "norm_"
      header = $1 OFS $2 OFS $3;
      for(i=4; i<=NF; i++){
         header = header OFS "norm_" $i;
      }
      print header;
      next;
    }
    {
      # Calculate the peak length (avoid division by zero)
      peak_length = $3 - $2;
      if (peak_length == 0) { peak_length = 1; }
      # Print the original BED coordinates
      printf("%s\t%s\t%s", $1, $2, $3);
      # For each BAM column (starting at field 4), compute normalized value
      for(i=4; i<=NF; i++){
        normalized = ($i * 1e9) / (peak_length * lib[i-3]);
        # Print with 3 decimal places
        printf("\t%.3f", normalized);
      }
      printf("\n");
    }' "$file" > "$out_file"
    
    echo "Normalized data saved to $out_file"
  fi
done
```
```R
library(tidyverse)

# Set working directory to where your *_multicov_normalized.txt files are stored
setwd("/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices_mm10/Fig_2D_R1_analyses/")

# List all normalized files
files <- list.files(pattern = "*_multicov_normalized.txt")

for(file in files){
  
  # Read the normalized file.
  # Assumes the first three columns are "chr", "start", "end"
  # and the remaining columns are normalized counts.
  dat <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Add a "peak" identifier for each row (each peak)
  dat <- dat %>% mutate(peak = row_number())
  
  # Identify sample columns: columns 4 through the last original column.
  sample_cols <- colnames(dat)[4:(ncol(dat)-1)]
  
  # Add 1 to all normalized values (to avoid zeros) in the sample columns.
  dat[, sample_cols] <- dat[, sample_cols] + 1
  
  # Pivot the sample columns to long format.
  # The column names are expected to be like:
  #   norm_ESC-MYC-G1, norm_ESC-YAP1-G1, ... , norm_EpiLC-d2-YAP1-G1
  dat_long <- dat %>% 
    pivot_longer(
      cols = sample_cols,
      names_to = "colName",
      values_to = "norm_val"
    ) %>%
    # Remove the "norm_" prefix
    mutate(colName = gsub("^norm_", "", colName)) %>%
    # Extract condition and TF factor.
    mutate(
      condition = case_when(
        grepl("^ESC", colName, ignore.case = TRUE)   ~ "ESC",
        grepl("^EpiLC", colName, ignore.case = TRUE) ~ "EpiLC",
        TRUE ~ NA_character_
      ),
      # Remove the condition tag (with separator . , - or _) and an optional "d2" for EpiLC.
      # Then convert to uppercase.
      factor = toupper(gsub("^(ESC[\\.-_]|EpiLC[\\.-_](d2[\\.-_])?)", "", colName, ignore.case = TRUE))
    ) %>%
    # Remove any trailing "-G1", ".G1", or "_G1" from the factor name
    mutate(factor = gsub("[-\\._]G1$", "", factor)) %>%
    filter(!is.na(condition))
  
  # For each peak and each TF, pivot so that we have one row with two columns: ESC and EpiLC.
  dat_wide <- dat_long %>%
    select(peak, factor, condition, norm_val) %>%
    pivot_wider(names_from = condition, values_from = norm_val)
  
  # Calculate the log2 ratio for each TF per peak: log2(ESC / EpiLC)
  dat_wide <- dat_wide %>%
    mutate(ratio = log2(ESC / EpiLC))
  
  # Reshape so that each row is a peak and each column is a TF.
  ratio_by_peak <- dat_wide %>%
    select(peak, factor, ratio) %>%
    pivot_wider(names_from = factor, values_from = ratio)
  
  # Reorder the columns to the desired order.
  desired_order <- c("OCT4", "SOX2", "NANOG", "MYC", "YAP1")
  ratio_by_peak <- ratio_by_peak %>% select(any_of(desired_order))
  
  # Write out the ratio file.
  out_file <- gsub("_multicov_normalized.txt", "_ratio.txt", file)
  write.table(ratio_by_peak, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("Processed file: ", file, " -> ratio file saved as: ", out_file)
}

logfile <- "session_info.log"
sink(logfile, append = TRUE)
sessionInfo()
if (requireNamespace("sessioninfo", quietly = TRUE)) sessioninfo::session_info()
sink()
message("Wrote session info to ", logfile)
```
```R
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
 [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       stringi_1.8.4     generics_0.1.3    glue_1.8.0        colorspace_2.1-1  hms_1.1.3        
[10] scales_1.3.0      grid_4.4.1        munsell_0.5.1     tzdb_0.4.0        lifecycle_1.0.4   compiler_4.4.1    sessioninfo_1.2.3 timechange_0.3.0  pkgconfig_2.0.3  
[19] rstudioapi_0.17.1 R6_2.5.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    tools_4.4.1       withr_3.0.2       gtable_0.3.6     
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.1 (2024-06-14)
 os       macOS Sonoma 14.7.2
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Berlin
 date     2025-05-02
 rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
 pandoc   NA
 quarto   1.4.555 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version date (UTC) lib source
 cli           3.6.3   2024-06-21 [1] CRAN (R 4.4.0)
 colorspace    2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
 crayon        1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
 dplyr       * 1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
 forcats     * 1.0.0   2023-01-29 [1] CRAN (R 4.4.0)
 generics      0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
 ggplot2     * 3.5.1   2024-04-23 [1] CRAN (R 4.4.0)
 glue          1.8.0   2024-09-30 [1] CRAN (R 4.4.1)
 gtable        0.3.6   2024-10-25 [1] CRAN (R 4.4.1)
 hms           1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
 lifecycle     1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
 lubridate   * 1.9.4   2024-12-08 [1] CRAN (R 4.4.1)
 magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
 munsell       0.5.1   2024-04-01 [1] CRAN (R 4.4.0)
 pillar        1.10.1  2025-01-07 [1] CRAN (R 4.4.1)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
 purrr       * 1.0.2   2023-08-10 [1] CRAN (R 4.4.0)
 R6            2.5.1   2021-08-19 [1] CRAN (R 4.4.0)
 readr       * 2.1.5   2024-01-10 [1] CRAN (R 4.4.0)
 rlang         1.1.4   2024-06-04 [1] CRAN (R 4.4.0)
 rstudioapi    0.17.1  2024-10-22 [1] CRAN (R 4.4.1)
 scales        1.3.0   2023-11-28 [1] CRAN (R 4.4.0)
 sessioninfo   1.2.3   2025-02-05 [1] CRAN (R 4.4.1)
 stringi       1.8.4   2024-05-06 [1] CRAN (R 4.4.0)
 stringr     * 1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
 tibble      * 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
 tidyr       * 1.3.1   2024-01-24 [1] CRAN (R 4.4.0)
 tidyselect    1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
 tidyverse   * 2.0.0   2023-02-22 [1] CRAN (R 4.4.0)
 timechange    0.3.0   2024-01-18 [1] CRAN (R 4.4.1)
 tzdb          0.4.0   2023-05-12 [1] CRAN (R 4.4.0)
 vctrs         0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
 withr         3.0.2   2024-10-28 [1] CRAN (R 4.4.1)

 [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
 * ── Packages attached to the search path.
```
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
