# DOA_ESC_TF_vs_TF
```bash
## Load Packages
library(edgeR)
library(ggplot2)
library(dplyr)

## Load Data and Assemble Count Matrices
rm(list = ls())

# Define the file path
file_path <- "/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices_mm10_ESC_only/"  # <-- CHANGE as needed

# ------------------------------------------------------------------------------
# 1) Reading ONLY ESC data for each of the 5 TFs: SOX2, OCT4, NANOG, MYC, YAP1
#    For each cell cycle phase (G1, S, G2), you have 2 replicates for each TF.
#    => 10 columns total per phase.
# ------------------------------------------------------------------------------
#
# We keep S1_G1 (for example) as a data frame with coordinate row names,
# and convert the other G1 replicates to named vectors. Then we cbind(...).

# ------------------------------
# -- G1 --
# ------------------------------
S1_G1 <- read.table(paste0(file_path, "ESC-SOX2-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)

S2_G1 <- read.table(paste0(file_path, "ESC-SOX2-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
S2_G1 <- S2_G1[, 1]

O1_G1 <- read.table(paste0(file_path, "ESC-OCT4-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O1_G1 <- O1_G1[, 1]

O2_G1 <- read.table(paste0(file_path, "ESC-OCT4-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O2_G1 <- O2_G1[, 1]

N1_G1 <- read.table(paste0(file_path, "ESC-NANOG-G1-1_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N1_G1 <- N1_G1[, 1]

N2_G1 <- read.table(paste0(file_path, "ESC-NANOG-G1-2_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N2_G1 <- N2_G1[, 1]

M1_G1 <- read.table(paste0(file_path, "ESC-MYC-G1-1_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M1_G1 <- M1_G1[, 1]

M2_G1 <- read.table(paste0(file_path, "ESC-MYC-G1-2_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M2_G1 <- M2_G1[, 1]

Y1_G1 <- read.table(paste0(file_path, "ESC-YAP1-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y1_G1 <- Y1_G1[, 1]

Y2_G1 <- read.table(paste0(file_path, "ESC-YAP1-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y2_G1 <- Y2_G1[, 1]

G1_merged <- cbind(
  S1_G1, 
  S2_G1, O1_G1, O2_G1, N1_G1, N2_G1, M1_G1, M2_G1, Y1_G1, Y2_G1
)
colnames(G1_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)

# ------------------------------
# -- S --
# ------------------------------
S1_S <- read.table(paste0(file_path, "ESC-SOX2-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)

S2_S <- read.table(paste0(file_path, "ESC-SOX2-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
S2_S <- S2_S[, 1]

O1_S <- read.table(paste0(file_path, "ESC-OCT4-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
O1_S <- O1_S[, 1]

O2_S <- read.table(paste0(file_path, "ESC-OCT4-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
O2_S <- O2_S[, 1]

N1_S <- read.table(paste0(file_path, "ESC-NANOG-S-1_counts.txt"), 
                   header=FALSE, row.names=1, check.names=FALSE)
N1_S <- N1_S[, 1]

N2_S <- read.table(paste0(file_path, "ESC-NANOG-S-2_counts.txt"), 
                   header=FALSE, row.names=1, check.names=FALSE)
N2_S <- N2_S[, 1]

M1_S <- read.table(paste0(file_path, "ESC-MYC-S-1_counts.txt"),   
                   header=FALSE, row.names=1, check.names=FALSE)
M1_S <- M1_S[, 1]

M2_S <- read.table(paste0(file_path, "ESC-MYC-S-2_counts.txt"),   
                   header=FALSE, row.names=1, check.names=FALSE)
M2_S <- M2_S[, 1]

Y1_S <- read.table(paste0(file_path, "ESC-YAP1-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
Y1_S <- Y1_S[, 1]

Y2_S <- read.table(paste0(file_path, "ESC-YAP1-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
Y2_S <- Y2_S[, 1]

S_merged <- cbind(
  S1_S,
  S2_S, O1_S, O2_S, N1_S, N2_S, M1_S, M2_S, Y1_S, Y2_S
)
colnames(S_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)

# ------------------------------
# -- G2 --
# ------------------------------
S1_G2 <- read.table(paste0(file_path, "ESC-SOX2-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)

S2_G2 <- read.table(paste0(file_path, "ESC-SOX2-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
S2_G2 <- S2_G2[, 1]

O1_G2 <- read.table(paste0(file_path, "ESC-OCT4-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O1_G2 <- O1_G2[, 1]

O2_G2 <- read.table(paste0(file_path, "ESC-OCT4-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O2_G2 <- O2_G2[, 1]

N1_G2 <- read.table(paste0(file_path, "ESC-NANOG-G2-1_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N1_G2 <- N1_G2[, 1]

N2_G2 <- read.table(paste0(file_path, "ESC-NANOG-G2-2_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N2_G2 <- N2_G2[, 1]

M1_G2 <- read.table(paste0(file_path, "ESC-MYC-G2-1_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M1_G2 <- M1_G2[, 1]

M2_G2 <- read.table(paste0(file_path, "ESC-MYC-G2-2_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M2_G2 <- M2_G2[, 1]

Y1_G2 <- read.table(paste0(file_path, "ESC-YAP1-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y1_G2 <- Y1_G2[, 1]

Y2_G2 <- read.table(paste0(file_path, "ESC-YAP1-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y2_G2 <- Y2_G2[, 1]

G2_merged <- cbind(
  S1_G2,
  S2_G2, O1_G2, O2_G2, N1_G2, N2_G2, M1_G2, M2_G2, Y1_G2, Y2_G2
)
colnames(G2_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)

## Setup edgeR Parameters
analyze_phase <- function(merged_data, phase, file_path) {
  
  # 5 groups: SOX2, OCT4, NANOG, MYC, YAP1
  group <- factor(c(
    rep("SOX2", 2),
    rep("OCT4", 2),
    rep("NANOG", 2),
    rep("MYC", 2),
    rep("YAP1", 2)
  ),
  levels = c("SOX2", "OCT4", "NANOG", "MYC", "YAP1"))
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # 10 pairwise contrasts among the 5 TFs
  # e.g. SOX2 - OCT4, SOX2 - NANOG, etc.
  contrast_matrix <- makeContrasts(
    SOX2_vs_OCT4   = SOX2 - OCT4,
    SOX2_vs_NANOG  = SOX2 - NANOG,
    SOX2_vs_MYC    = SOX2 - MYC,
    SOX2_vs_YAP1   = SOX2 - YAP1,
    OCT4_vs_NANOG  = OCT4 - NANOG,
    OCT4_vs_MYC    = OCT4 - MYC,
    OCT4_vs_YAP1   = OCT4 - YAP1,
    NANOG_vs_MYC   = NANOG - MYC,
    NANOG_vs_YAP1  = NANOG - YAP1,
    MYC_vs_YAP1    = MYC - YAP1,
    levels = design
  )
  
  # edgeR Differential Binding Analysis
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Plot MDS
  pdf(file = paste0(file_path, "MDS_Plot_", phase, ".pdf"))
  plotMDS(dge, main = paste("MDS Plot -", phase), col = as.numeric(group))
  dev.off()
  
  # Fit the model
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # We'll loop over each contrast to generate separate results, plots, and BEDs
  contrast_names <- colnames(contrast_matrix)
  
  for (cn in contrast_names) {
    
    # Perform LRT for this specific contrast
    lrt <- glmLRT(fit, contrast = contrast_matrix[, cn])
    results <- topTags(lrt, n = Inf)$table
    results$FDR <- p.adjust(results$PValue, method = "BH")
    
    # Label up/down/not significant
    sig_regions <- results %>%
      mutate(
        Category = case_when(
          logFC < -0.5 & FDR < 0.05 ~ "Downregulated",
          logFC > 0.5 & FDR < 0.05 ~ "Upregulated",
          TRUE ~ "Non-significant"
        )
      )
    
    # Volcano Plot
    category_counts <- sig_regions %>%
      group_by(Category) %>%
      summarize(Count = n(), .groups = "drop")
    legend_labels <- category_counts %>%
      mutate(Label = paste(Category, "(", Count, ")")) %>%
      pull(Label)
    
    volcano_plot <- ggplot(data = sig_regions, aes(x = logFC, y = -log10(FDR), color = Category)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c("Downregulated" = "darkgrey", 
                                    "Non-significant" = "lightblue", 
                                    "Upregulated" = "lightgrey"),
                         labels = legend_labels) +
      theme_classic() +
      ggtitle(paste0(cn, " - ", phase, " Phase")) +
      theme(legend.position = "right") +
      xlim(c(-max(abs(results$logFC)), max(abs(results$logFC)))) +
      labs(x = "Log Fold Change", y = "-Log10 FDR", color = "Differential Expression")
    
    ggsave(filename = paste0(file_path, cn, "_Volcano_Plot_", phase, ".pdf"), plot = volcano_plot)
    
    # BED output (Up & Down)
    results$chr   <- sapply(strsplit(rownames(results), "[:-]"), `[`, 1)
    results$start <- as.numeric(sapply(strsplit(rownames(results), "[:-]"), `[`, 2))
    results$end   <- as.numeric(sapply(strsplit(rownames(results), "[:-]"), `[`, 3))
    
    down_regions <- results %>%
      filter(logFC < -0.5 & FDR < 0.05) %>%
      select(chr, start, end) %>%
      mutate(category = "DOWN")
    
    up_regions <- results %>%
      filter(logFC > 0.5 & FDR < 0.05) %>%
      select(chr, start, end) %>%
      mutate(category = "UP")
    
    combined_regions <- bind_rows(down_regions, up_regions)
    bed_file <- paste0(file_path, cn, "_edgeR_DBRs_", phase, ".bed")
    write.table(combined_regions, file = bed_file, quote = FALSE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
    
    # (Optional) If you want separate UP & DOWN bed files, you can
    # replicate the approach from before. For example:
    #
    up_bed_file <- paste0(file_path, cn, "_edgeR_UP_", phase, ".bed")
    write.table(up_regions, file=up_bed_file, quote=FALSE, sep="\t", 
                 row.names=FALSE, col.names=FALSE)
  
    down_bed_file <- paste0(file_path, cn, "_edgeR_DOWN_", phase, ".bed")
    write.table(down_regions, file=down_bed_file, quote=FALSE, sep="\t",
                 row.names=FALSE, col.names=FALSE)
    
    # Bar Plot: Lost vs Gained
    lost <- nrow(down_regions)
    gained <- nrow(up_regions)
    count_data <- data.frame(
      Category = c("Lost Regions", "Gained Regions"),
      Count = c(lost, gained)
    ) %>%
      mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))
    
    region_count_plot <- ggplot(count_data, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
      labs(title = paste("Number of Regions:", cn, phase),
           x = "", y = "Number of Regions", fill = "Differential Binding") +
      theme_minimal() +
      scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))
    
    ggsave(filename = paste0(file_path, cn, "_Regions_", phase, ".pdf"), plot = region_count_plot)
    
  } # end of loop over contrasts
}

# ----------
# 3) Run the analysis for each phase (G1, S, G2)
# ----------
analyze_phase(G1_merged, "G1", file_path)
analyze_phase(S_merged, "S", file_path)
analyze_phase(G2_merged, "G2", file_path)

# ----------
# 4) (Optional) Summary Table
#    This code tallies how many Down / Up for each contrast & phase.
# ----------
analyze_phase_summary.table <- function(merged_data) {
  
  group <- factor(c(
    rep("SOX2", 2),
    rep("OCT4", 2),
    rep("NANOG", 2),
    rep("MYC", 2),
    rep("YAP1", 2)
  ),
  levels = c("SOX2", "OCT4", "NANOG", "MYC", "YAP1"))
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # Same 10 pairwise contrasts
  contr <- makeContrasts(
    SOX2_vs_OCT4   = SOX2 - OCT4,
    SOX2_vs_NANOG  = SOX2 - NANOG,
    SOX2_vs_MYC    = SOX2 - MYC,
    SOX2_vs_YAP1   = SOX2 - YAP1,
    OCT4_vs_NANOG  = OCT4 - NANOG,
    OCT4_vs_MYC    = OCT4 - MYC,
    OCT4_vs_YAP1   = OCT4 - YAP1,
    NANOG_vs_MYC   = NANOG - MYC,
    NANOG_vs_YAP1  = NANOG - YAP1,
    MYC_vs_YAP1    = MYC - YAP1,
    levels = design
  )
  
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # For each of the 10 contrasts, compute up/down
  out_list <- list()
  for (cn in colnames(contr)) {
    lrt <- glmLRT(fit, contrast=contr[, cn])
    results <- topTags(lrt, n = Inf)$table
    results$FDR <- p.adjust(results$PValue, method = "BH")
    
    down_count <- sum(results$logFC < -0.5 & results$FDR < 0.05)
    up_count   <- sum(results$logFC >  0.5 & results$FDR < 0.05)
    
    out_list[[cn]] <- data.frame(
      Contrast = cn,
      Lost     = down_count,
      Gained   = up_count
    )
  }
  do.call(rbind, out_list)
}

# Build final summary over phases
all_summary <- data.frame()
for (ph in c("G1", "S", "G2")) {
  merged_mat <- switch(ph,
                       "G1" = G1_merged,
                       "S"  = S_merged,
                       "G2" = G2_merged)
  res <- analyze_phase_summary.table(merged_mat)
  res$Phase <- ph
  all_summary <- rbind(all_summary, res)
}

# Reorder columns to something like: Phase, Contrast, Lost, Gained
all_summary <- all_summary[, c("Phase", "Contrast", "Lost", "Gained")]

# Write out summary table
output_summary_file <- paste0(file_path, "ESC_All5TFs_Summary_Table.csv")
write.table(all_summary, file = output_summary_file, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

message("Summary table saved to: ", output_summary_file)

logfile <- "session_info.log"
sink(logfile, append = TRUE)
sessionInfo()
if (requireNamespace("sessioninfo", quietly = TRUE)) sessioninfo::session_info()
sink()
message("Wrote session info to ", logfile)
```
# Session info DOA_ESC_TF_vs_TF
```bash
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

loaded via a namespace (and not attached):
[1] compiler_4.4.1    tools_4.4.1       rstudioapi_0.17.1
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
 [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       generics_0.1.3    textshaping_1.0.0 labeling_0.4.3    glue_1.8.0        statmod_1.5.0    
 [9] colorspace_2.1-1  ragg_1.3.3        locfit_1.5-9.10   scales_1.3.0      grid_4.4.1        munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4  
[17] compiler_4.4.1    Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2      systemfonts_1.2.0 lattice_0.22-6    R6_2.5.1         
[25] tidyselect_1.2.1  pillar_1.10.1     splines_4.4.1     magrittr_2.0.3    tools_4.4.1       withr_3.0.2       gtable_0.3.6
```
# DOA_ESC_SNO_vs_TF
```bash
## Load Packages
library(edgeR)
library(ggplot2)
library(dplyr)

## Load Data and Assemble Count Matrices
rm(list = ls())

# Define the file path
file_path <- "/Users/hansel01/Desktop/Desktop_2/job_application_082016/CMMC/CMMC_RHH.lab/CMMC_Projects/DynaTag/seq_data_DynaTag/DynaTag/ESC_EpiLC_DynaTag/TF_count_matrices_mm10_ESC_only/"  

# ------------------------------------------------------------------------------
# 1) Reading ONLY ESC data for each TF: SOX2, OCT4, NANOG, MYC, YAP1
#    For each cell cycle phase (G1, S, G2), you have 2 replicates for each TF.
#    SNO = SOX2 + OCT4 + NANOG => total 6 columns
#    MYC => 2 columns
#    YAP1 => 2 columns
#    => 10 columns total per phase.
# ------------------------------------------------------------------------------
# We keep one file (S1_G1, etc.) as a data frame to preserve row names, and
# convert all other files in that phase to named vectors. Then cbind() them.

# ------------------------------
# -- G1 --
# ------------------------------
S1_G1 <- read.table(paste0(file_path, "ESC-SOX2-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)

S2_G1 <- read.table(paste0(file_path, "ESC-SOX2-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
S2_G1 <- S2_G1[, 1]

O1_G1 <- read.table(paste0(file_path, "ESC-OCT4-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O1_G1 <- O1_G1[, 1]

O2_G1 <- read.table(paste0(file_path, "ESC-OCT4-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O2_G1 <- O2_G1[, 1]

N1_G1 <- read.table(paste0(file_path, "ESC-NANOG-G1-1_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N1_G1 <- N1_G1[, 1]

N2_G1 <- read.table(paste0(file_path, "ESC-NANOG-G1-2_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N2_G1 <- N2_G1[, 1]

M1_G1 <- read.table(paste0(file_path, "ESC-MYC-G1-1_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M1_G1 <- M1_G1[, 1]

M2_G1 <- read.table(paste0(file_path, "ESC-MYC-G1-2_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M2_G1 <- M2_G1[, 1]

Y1_G1 <- read.table(paste0(file_path, "ESC-YAP1-G1-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y1_G1 <- Y1_G1[, 1]

Y2_G1 <- read.table(paste0(file_path, "ESC-YAP1-G1-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y2_G1 <- Y2_G1[, 1]

G1_merged <- cbind(
  S1_G1, 
  S2_G1, O1_G1, O2_G1, N1_G1, N2_G1, M1_G1, M2_G1, Y1_G1, Y2_G1
)
colnames(G1_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)

# ------------------------------
# -- S --
# ------------------------------
S1_S <- read.table(paste0(file_path, "ESC-SOX2-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)

S2_S <- read.table(paste0(file_path, "ESC-SOX2-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
S2_S <- S2_S[, 1]

O1_S <- read.table(paste0(file_path, "ESC-OCT4-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
O1_S <- O1_S[, 1]

O2_S <- read.table(paste0(file_path, "ESC-OCT4-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
O2_S <- O2_S[, 1]

N1_S <- read.table(paste0(file_path, "ESC-NANOG-S-1_counts.txt"), 
                   header=FALSE, row.names=1, check.names=FALSE)
N1_S <- N1_S[, 1]

N2_S <- read.table(paste0(file_path, "ESC-NANOG-S-2_counts.txt"), 
                   header=FALSE, row.names=1, check.names=FALSE)
N2_S <- N2_S[, 1]

M1_S <- read.table(paste0(file_path, "ESC-MYC-S-1_counts.txt"),   
                   header=FALSE, row.names=1, check.names=FALSE)
M1_S <- M1_S[, 1]

M2_S <- read.table(paste0(file_path, "ESC-MYC-S-2_counts.txt"),   
                   header=FALSE, row.names=1, check.names=FALSE)
M2_S <- M2_S[, 1]

Y1_S <- read.table(paste0(file_path, "ESC-YAP1-S-1_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
Y1_S <- Y1_S[, 1]

Y2_S <- read.table(paste0(file_path, "ESC-YAP1-S-2_counts.txt"),  
                   header=FALSE, row.names=1, check.names=FALSE)
Y2_S <- Y2_S[, 1]

S_merged <- cbind(
  S1_S,
  S2_S, O1_S, O2_S, N1_S, N2_S, M1_S, M2_S, Y1_S, Y2_S
)
colnames(S_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)

# ------------------------------
# -- G2 --
# ------------------------------
S1_G2 <- read.table(paste0(file_path, "ESC-SOX2-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)

S2_G2 <- read.table(paste0(file_path, "ESC-SOX2-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
S2_G2 <- S2_G2[, 1]

O1_G2 <- read.table(paste0(file_path, "ESC-OCT4-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O1_G2 <- O1_G2[, 1]

O2_G2 <- read.table(paste0(file_path, "ESC-OCT4-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
O2_G2 <- O2_G2[, 1]

N1_G2 <- read.table(paste0(file_path, "ESC-NANOG-G2-1_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N1_G2 <- N1_G2[, 1]

N2_G2 <- read.table(paste0(file_path, "ESC-NANOG-G2-2_counts.txt"), 
                    header=FALSE, row.names=1, check.names=FALSE)
N2_G2 <- N2_G2[, 1]

M1_G2 <- read.table(paste0(file_path, "ESC-MYC-G2-1_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M1_G2 <- M1_G2[, 1]

M2_G2 <- read.table(paste0(file_path, "ESC-MYC-G2-2_counts.txt"),   
                    header=FALSE, row.names=1, check.names=FALSE)
M2_G2 <- M2_G2[, 1]

Y1_G2 <- read.table(paste0(file_path, "ESC-YAP1-G2-1_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y1_G2 <- Y1_G2[, 1]

Y2_G2 <- read.table(paste0(file_path, "ESC-YAP1-G2-2_counts.txt"),  
                    header=FALSE, row.names=1, check.names=FALSE)
Y2_G2 <- Y2_G2[, 1]

G2_merged <- cbind(
  S1_G2,
  S2_G2, O1_G2, O2_G2, N1_G2, N2_G2, M1_G2, M2_G2, Y1_G2, Y2_G2
)
colnames(G2_merged) <- c(
  "SOX2_1","SOX2_2",
  "OCT4_1","OCT4_2",
  "NANOG_1","NANOG_2",
  "MYC_1","MYC_2",
  "YAP1_1","YAP1_2"
)


## Setup edgeR Parameters
analyze_phase <- function(merged_data, phase, file_path) {
  # ------
  # Create a multi-level factor: 6 "SNO" + 2 "MYC" + 2 "YAP1"
  # We'll make 3 contrasts:
  #  1) SNO vs YAP1
  #  2) SNO vs MYC
  #  3) MYC vs YAP1
  # ------
  
  group <- factor(c(
    rep("SNO", 6), 
    rep("MYC", 2), 
    rep("YAP1", 2)
  ),
  levels = c("SNO", "MYC", "YAP1"))
  
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  contrast_matrix <- makeContrasts(
    SNO_vs_MYC   = SNO - MYC,
    SNO_vs_YAP1  = SNO - YAP1,
    MYC_vs_YAP1  = MYC - YAP1,
    levels = design
  )
  
  # edgeR Differential Binding Analysis
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Plot MDS
  pdf(file = paste0(file_path, "MDS_Plot_", phase, ".pdf"))
  plotMDS(dge, main = paste("MDS Plot -", phase), col = as.numeric(group))
  dev.off()
  
  # Fit the model
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # We'll loop over each contrast to generate separate results, plots, and BEDs
  contrast_names <- colnames(contrast_matrix)
  
  for (cn in contrast_names) {
    
    # Perform LRT for this specific contrast
    lrt <- glmLRT(fit, contrast = contrast_matrix[, cn])
    results <- topTags(lrt, n = Inf)$table
    results$FDR <- p.adjust(results$PValue, method = "BH")
    
    # Label up/down/not significant
    sig_regions <- results %>%
      mutate(
        Category = case_when(
          logFC < -0.5 & FDR < 0.05 ~ "Downregulated",
          logFC > 0.5 & FDR < 0.05 ~ "Upregulated",
          TRUE ~ "Non-significant"
        )
      )
    
    # Volcano Plot
    category_counts <- sig_regions %>%
      group_by(Category) %>%
      summarize(Count = n(), .groups = "drop")
    legend_labels <- category_counts %>%
      mutate(Label = paste(Category, "(", Count, ")")) %>%
      pull(Label)
    
    volcano_plot <- ggplot(data = sig_regions, aes(x = logFC, y = -log10(FDR), color = Category)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c("Downregulated" = "darkgrey", 
                                    "Non-significant" = "lightblue", 
                                    "Upregulated" = "lightgrey"),
                         labels = legend_labels) +
      theme_classic() +
      ggtitle(paste0(cn, " - ", phase, " Phase")) +
      theme(legend.position = "right") +
      xlim(c(-max(abs(results$logFC)), max(abs(results$logFC)))) +
      labs(x = "Log Fold Change", y = "-Log10 FDR", color = "Differential Expression")
    
    ggsave(filename = paste0(file_path, cn, "_Volcano_Plot_", phase, ".pdf"), plot = volcano_plot)
    
    # ## Updated BED file output: Combined, UP, and DOWN ##
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
    
    # 1) Combined: UP + DOWN
    combined_regions <- bind_rows(down_regions, up_regions)
    bed_file <- paste0(file_path, cn, "_edgeR_DBRs_", phase, ".bed")
    write.table(combined_regions, file = bed_file, quote = FALSE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
    
    # 2) UP-only
    up_bed_file <- paste0(file_path, cn, "_edgeR_UP_", phase, ".bed")
    write.table(up_regions, file = up_bed_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    # 3) DOWN-only
    down_bed_file <- paste0(file_path, cn, "_edgeR_DOWN_", phase, ".bed")
    write.table(down_regions, file = down_bed_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    # Bar Plot: Lost vs Gained
    lost <- nrow(down_regions)
    gained <- nrow(up_regions)
    count_data <- data.frame(
      Category = c("Lost Regions", "Gained Regions"),
      Count = c(lost, gained)
    ) %>%
      mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))
    
    region_count_plot <- ggplot(count_data, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
      labs(title = paste("Number of Regions:", cn, phase),
           x = "", y = "Number of Regions", fill = "Differential Binding") +
      theme_minimal() +
      scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))
    
    ggsave(filename = paste0(file_path, cn, "_Regions_", phase, ".pdf"), plot = region_count_plot)
    
  } # end of loop over contrasts
}


# ----------
# 3) Run the analysis for each phase (G1, S, G2)
# ----------
analyze_phase(G1_merged, "G1", file_path)
analyze_phase(S_merged, "S", file_path)
analyze_phase(G2_merged, "G2", file_path)

# ----------
# 4) (Optional) Summary Table
# ----------
analyze_phase_summary.table <- function(merged_data) {
  group <- factor(c(rep("SNO", 6), rep("MYC", 2), rep("YAP1", 2)),
                  levels = c("SNO", "MYC", "YAP1"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  contr <- makeContrasts(
    SNO_vs_MYC   = SNO - MYC,
    SNO_vs_YAP1  = SNO - YAP1,
    MYC_vs_YAP1  = MYC - YAP1,
    levels = design
  )
  dge <- DGEList(counts = merged_data)
  dge$samples$group <- group
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # For each of the three contrasts, compute up/down
  out_list <- list()
  for (cn in colnames(contr)) {
    lrt <- glmLRT(fit, contrast=contr[, cn])
    results <- topTags(lrt, n = Inf)$table
    results$FDR <- p.adjust(results$PValue, method = "BH")
    down_count <- sum(results$logFC < -0.5 & results$FDR < 0.05)
    up_count   <- sum(results$logFC >  0.5 & results$FDR < 0.05)
    
    out_list[[cn]] <- data.frame(
      Contrast = cn,
      Lost     = down_count,
      Gained   = up_count
    )
  }
  do.call(rbind, out_list)
}

# Build final summary over phases
all_summary <- data.frame()
for (ph in c("G1", "S", "G2")) {
  merged_mat <- switch(ph,
                       "G1" = G1_merged,
                       "S"  = S_merged,
                       "G2" = G2_merged)
  res <- analyze_phase_summary.table(merged_mat)
  res$Phase <- ph
  all_summary <- rbind(all_summary, res)
}

all_summary <- all_summary[, c("Phase", "Contrast", "Lost", "Gained")]

output_summary_file <- paste0(file_path, "ESC_SNO_MYC_YAP1_Summary_Table.csv")
write.table(all_summary, file = output_summary_file, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
message("Summary table saved to: ", output_summary_file)

logfile <- "session_info.log"
sink(logfile, append = TRUE)
sessionInfo()
if (requireNamespace("sessioninfo", quietly = TRUE)) sessioninfo::session_info()
sink()
message("Wrote session info to ", logfile)
```
# Session info DOA_ESC_SNO_vs_TF
```bash
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

loaded via a namespace (and not attached):
[1] compiler_4.4.1    tools_4.4.1       rstudioapi_0.17.1
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
 [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       generics_0.1.3    textshaping_1.0.0 labeling_0.4.3    glue_1.8.0        statmod_1.5.0    
 [9] colorspace_2.1-1  ragg_1.3.3        locfit_1.5-9.10   scales_1.3.0      grid_4.4.1        munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4  
[17] compiler_4.4.1    Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2      systemfonts_1.2.0 lattice_0.22-6    R6_2.5.1         
[25] tidyselect_1.2.1  pillar_1.10.1     splines_4.4.1     magrittr_2.0.3    tools_4.4.1       withr_3.0.2       gtable_0.3.6     
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
 [1] vctrs_0.6.5       cli_3.6.3         rlang_1.1.4       generics_0.1.3    textshaping_1.0.0 labeling_0.4.3    glue_1.8.0        statmod_1.5.0    
 [9] colorspace_2.1-1  ragg_1.3.3        locfit_1.5-9.10   scales_1.3.0      grid_4.4.1        munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4  
[17] compiler_4.4.1    sessioninfo_1.2.3 Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2      systemfonts_1.2.0 lattice_0.22-6   
[25] R6_2.5.1          tidyselect_1.2.1  pillar_1.10.1     splines_4.4.1     magrittr_2.0.3    tools_4.4.1       withr_3.0.2       gtable_0.3.6     
─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version  date (UTC) lib source
 cli           3.6.3    2024-06-21 [1] CRAN (R 4.4.0)
 colorspace    2.1-1    2024-07-26 [1] CRAN (R 4.4.0)
 dplyr       * 1.1.4    2023-11-17 [1] CRAN (R 4.4.0)
 edgeR       * 4.2.2    2024-10-13 [1] Bioconductor 3.19 (R 4.4.1)
 farver        2.1.2    2024-05-13 [1] CRAN (R 4.4.0)
 generics      0.1.3    2022-07-05 [1] CRAN (R 4.4.0)
 ggplot2     * 3.5.1    2024-04-23 [1] CRAN (R 4.4.0)
 glue          1.8.0    2024-09-30 [1] CRAN (R 4.4.1)
 gtable        0.3.6    2024-10-25 [1] CRAN (R 4.4.1)
 labeling      0.4.3    2023-08-29 [1] CRAN (R 4.4.0)
 lattice       0.22-6   2024-03-20 [1] CRAN (R 4.4.1)
 lifecycle     1.0.4    2023-11-07 [1] CRAN (R 4.4.0)
 limma       * 3.60.6   2024-10-02 [1] Bioconductor 3.19 (R 4.4.1)
 locfit        1.5-9.10 2024-06-24 [1] CRAN (R 4.4.0)
 magrittr      2.0.3    2022-03-30 [1] CRAN (R 4.4.0)
 munsell       0.5.1    2024-04-01 [1] CRAN (R 4.4.0)
 pillar        1.10.1   2025-01-07 [1] CRAN (R 4.4.1)
 pkgconfig     2.0.3    2019-09-22 [1] CRAN (R 4.4.0)
 R6            2.5.1    2021-08-19 [1] CRAN (R 4.4.0)
 ragg          1.3.3    2024-09-11 [1] CRAN (R 4.4.1)
 Rcpp          1.0.14   2025-01-12 [1] CRAN (R 4.4.1)
 rlang         1.1.4    2024-06-04 [1] CRAN (R 4.4.0)
 rstudioapi    0.17.1   2024-10-22 [1] CRAN (R 4.4.1)
 scales        1.3.0    2023-11-28 [1] CRAN (R 4.4.0)
 sessioninfo   1.2.3    2025-02-05 [1] CRAN (R 4.4.1)
 statmod       1.5.0    2023-01-06 [1] CRAN (R 4.4.0)
 systemfonts   1.2.0    2025-01-16 [1] CRAN (R 4.4.1)
 textshaping   1.0.0    2025-01-20 [1] CRAN (R 4.4.1)
 tibble        3.2.1    2023-03-20 [1] CRAN (R 4.4.0)
 tidyselect    1.2.1    2024-03-11 [1] CRAN (R 4.4.0)
 vctrs         0.6.5    2023-12-01 [1] CRAN (R 4.4.0)
 withr         3.0.2    2024-10-28 [1] CRAN (R 4.4.1)

 [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
 * ── Packages attached to the search path.

───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

```
