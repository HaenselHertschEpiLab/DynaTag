# Preprocessing in DynaTag_sc_ESC.md
# Plot FRiP Scores as Violin in R
```R
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set the directory path
frip_path <- "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/snDynaTag_ICELL8/FRiP/"

# Read data from each file
OCT4 <- read.table(paste0(frip_path, "FRiP_scores_OCT4.txt"))$V6
NANOG <- read.table(paste0(frip_path, "FRiP_scores_NANOG.txt"))$V6
MYC <- read.table(paste0(frip_path, "FRiP_scores_MYC.txt"))$V6
YAP1 <- read.table(paste0(frip_path, "FRiP_scores_YAP1.txt"))$V6

# Create data frames with a common column for grouping
OCT4_df <- data.frame(Sample = rep("OCT4", length(OCT4)), Score = OCT4)
NANOG_df <- data.frame(Sample = rep("NANOG", length(NANOG)), Score = NANOG)
MYC_df <- data.frame(Sample = rep("MYC", length(MYC)), Score = MYC)
YAP1_df <- data.frame(Sample = rep("YAP1", length(YAP1)), Score = YAP1)

# Bind data frames together
combined_data <- bind_rows(OCT4_df, NANOG_df, MYC_df, YAP1_df)

# Create violin plot with facets, customized y-axis limits, and no x-axis label
ggplot(combined_data, aes(x = "", y = Score)) +
  geom_violin(fill = "lightblue", color = "blue") +
  labs(title = "FRiP Scores Violin Plot", x = NULL, y = "FRiP Score") +
  facet_wrap(~ Sample, scales = "free") +
  ylim(0, 1) +  # Set y-axis limits from 0 to 1
  theme_minimal() +
  xlab(NULL)  # Remove x-axis label

```
