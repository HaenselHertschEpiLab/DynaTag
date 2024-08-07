# Preprocessing in DynaTag_bulk_ESC_EpiLC.md
# Differential Binding Analysis
## Load Packages
```R
library(edgeR)
library(ggplot2)

```
## Load Data and Assemble Count Matrices
```R
rm(list = ls())

# Define the file path
file_path <- "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/count_matrices/"

# Define the TF to analyze
TF <- "OCT4"  # Change this to any other TF you want to analyze

# Read count data from files
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

# Check if row names are identical in all matrices for G1 phase
identical_row_names_G1 <- all(
  identical(rownames(ESC_G1_1), rownames(ESC_G1_2)),
  identical(rownames(ESC_G1_1), rownames(EpiLC_G1_1)),
  identical(rownames(ESC_G1_1), rownames(EpiLC_G1_2))
)

# Check if row names are identical in all matrices for S phase
identical_row_names_S <- all(
  identical(rownames(ESC_S_1), rownames(ESC_S_2)),
  identical(rownames(ESC_S_1), rownames(EpiLC_S_1)),
  identical(rownames(ESC_S_1), rownames(EpiLC_S_2))
)

# Check if row names are identical in all matrices for G2 phase
identical_row_names_G2 <- all(
  identical(rownames(ESC_G2_1), rownames(ESC_G2_2)),
  identical(rownames(ESC_G2_1), rownames(EpiLC_G2_1)),
  identical(rownames(ESC_G2_1), rownames(EpiLC_G2_2))
)


ESC_G1_1 <- ESC_G1_1
ESC_G1_2 <- ESC_G1_2[, 1]
EpiLC_G1_1 <- EpiLC_G1_1[, 1]
EpiLC_G1_2 <- EpiLC_G1_2[, 1]

ESC_S_1 <- ESC_S_1
ESC_S_2 <- ESC_S_2[, 1]
EpiLC_S_1 <- EpiLC_S_1[, 1]
EpiLC_S_2 <- EpiLC_S_2[, 1]

ESC_G2_1 <- ESC_G2_1
ESC_G2_2 <- ESC_G2_2[, 1]
EpiLC_G2_1 <- EpiLC_G2_1[, 1]
EpiLC_G2_2 <- EpiLC_G2_2[, 1]


# Merging Data
G1_merged <- cbind(ESC_G1_1, ESC_G1_2, EpiLC_G1_1, EpiLC_G1_2)
colnames(G1_merged)[1] <- "ESC_G1_1"

S_merged <- cbind(ESC_S_1, ESC_S_2, EpiLC_S_1, EpiLC_S_2)
colnames(S_merged)[1] <- "ESC_S_1"

G2_merged <- cbind(ESC_G2_1, ESC_G2_2, EpiLC_G2_1, EpiLC_G2_2)
colnames(G2_merged)[1] <- "ESC_G2_1"
```
## Setup edgeR Parametres
```R
# Set the reference group manually (change "PBS" to the group you want as the reference)
reference_group <- "ESC"

# Create a factor variable for group and specify the reference level
group <- factor(c("ESC", "ESC", "EpiLC", "EpiLC"), levels = c("ESC", "EpiLC"))

# Design Matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Create contrast with ESC as the reference group
contrast <- makeContrasts(reference_vs_other = paste0("EpiLC - ", reference_group), levels = design)
```
## Differential Binding Analysis per Cell Cycle Phase
```R
# Data Preprocessing
dge_G1 <- DGEList(counts = G1_merged)
dge_G1$samples$group <- group  # Adding group information to DGEList object
dge_G1 <- calcNormFactors(dge_G1, method = "TMM")

# Plot MDS
plotMDS(dge_G1, main="", col=as.numeric(group)) # Coloring by group

dge_G1 <- estimateGLMCommonDisp(dge_G1, design)
dge_G1 <- estimateGLMTagwiseDisp(dge_G1, design)

# Estimate Dispersion
dge_G1 <- estimateDisp(dge_G1)

# Fit GLM
fit_G1 <- glmFit(dge_G1, design)

# Likelihood Ratio Test (LRT)
lrt_G1 <- glmLRT(fit_G1, contrast = contrast)

# Extract and adjust results
results_G1 <- topTags(lrt_G1, n = Inf)
results_G1 <- results_G1$table
results_G1$FDR <- p.adjust(results_G1$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_G1)

# Volcano Plot
res_G1 <- na.omit(results_G1)
sig_G1 <- sum(res_G1$FDR < 0.05)
UP_G1 <- sum(res_G1$logFC > 0.5 & res_G1$FDR < 0.05)
DOWN_G1 <- sum(res_G1$logFC < -0.5 & res_G1$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_G1, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_G1, ")"))
}
ggplot(data = res_G1, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_G1$logFC)), max(abs(res_G1$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "G0/G1 Phase"))




# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_G1 <- DOWN_G1

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_G1 <- UP_G1

# Create a data frame for the counts
count_data_G1 <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_G1, gained_regions_G1)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_G1 <- count_data_G1 %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_G1 <- ggplot(count_data_G1, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in G0/G1 Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_G1)

ggsave(filename = paste0(file_path, TF, "_Regions_G1.pdf"), plot = plot_G1)



###### S #####

# Data Preprocessing
dge_S <- DGEList(counts = S_merged)
dge_S$samples$group <- group  # Adding group information to DGEList object
dge_S <- calcNormFactors(dge_S, method = "TMM")

# Plot MDS
plotMDS(dge_S, main="", col=as.numeric(group)) # Coloring by group

dge_S <- estimateGLMCommonDisp(dge_S, design)
dge_S <- estimateGLMTagwiseDisp(dge_S, design)

# Estimate Dispersion
dge_S <- estimateDisp(dge_S)

# Fit GLM
fit_S <- glmFit(dge_S, design)

# Likelihood Ratio Test (LRT)
lrt_S <- glmLRT(fit_S, contrast = contrast)

# Extract and adjust results
results_S <- topTags(lrt_S, n = Inf)
results_S <- results_S$table
results_S$FDR <- p.adjust(results_S$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_S)

# Volcano Plot
res_S <- na.omit(results_S)
sig_S <- sum(res_S$FDR < 0.05)
UP_S <- sum(res_S$logFC > 0.5 & res_S$FDR < 0.05)
DOWN_S <- sum(res_S$logFC < -0.5 & res_S$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_S, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_S, ")"))
}
ggplot(data = res_S, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_S$logFC)), max(abs(res_S$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "S Phase"))
###

# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_S <- DOWN_S

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_S <- UP_S

# Create a data frame for the counts
count_data_S <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_S, gained_regions_S)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_S <- count_data_S %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_S <- ggplot(count_data_S, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in S Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_S)

ggsave(filename = paste0(file_path, TF, "_Regions_S.pdf"), plot = plot_S)



###### G2 #####

# Data Preprocessing
dge_G2 <- DGEList(counts = G2_merged)
dge_G2$samples$group <- group  # Adding group information to DGEList object
dge_G2 <- calcNormFactors(dge_G2, method = "TMM")

# Plot MDS
plotMDS(dge_G2, main="", col=as.numeric(group)) # Coloring by group

dge_G2 <- estimateGLMCommonDisp(dge_G2, design)
dge_G2 <- estimateGLMTagwiseDisp(dge_G2, design)

# Estimate Dispersion
dge_G2 <- estimateDisp(dge_G2)

# Fit GLM
fit_G2 <- glmFit(dge_G2, design)

# Likelihood Ratio Test (LRT)
lrt_G2 <- glmLRT(fit_G2, contrast = contrast)

# Extract and adjust results
results_G2 <- topTags(lrt_G2, n = Inf)
results_G2 <- results_G2$table
results_G2$FDR <- p.adjust(results_G2$PValue, method = "BH") # Adjust p-values using the Benjamini-Hochberg method

# Plot Biological Coefficient of Variation (BCV)
plotBCV(dge_G2)

#Volcano Plot
res_G2 <- na.omit(results_G2)
sig_G2 <- sum(res_G2$FDR < 0.05)
UP_G2 <- sum(res_G2$logFC > 0.5 & res_G2$FDR < 0.05)
DOWN_G2 <- sum(res_G2$logFC < -0.5 & res_G2$FDR < 0.05)
if (sig_G1 == 0) {
  legend_labels <- c("Non-significant")
} else {
  legend_labels <- c(paste("Downregulated (", DOWN_G2, ")"),
                     "Non-significant",
                     paste("Upregulated (", UP_G2, ")"))
}
ggplot(data = res_G2, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(col = ifelse(logFC < -0.5 & FDR < 0.05, "darkgrey",
                              ifelse(logFC > 0.5 & FDR < 0.05, "lightgrey", "lightblue"))),
             alpha = 0.5) +
  scale_color_manual(values = c("darkgrey", "lightgrey", "lightblue"),
                     labels = legend_labels,
                     name = "Differential Expression") +
  theme_classic() +
  xlim(c(-max(abs(res_G2$logFC)), max(abs(res_G2$logFC)))) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20),
        panel.border = element_rect(color = "black", fill = NA)) +
  ggtitle(paste(TF, "G2/M Phase"))
###

# Count the number of regions where logFC is less than 0 (indicating loss in Chemo)
lost_regions_G2 <- DOWN_G2

# Count the number of regions where logFC is greater than 0 (indicating gain in Chemo)
gained_regions_G2 <- UP_G2

# Create a data frame for the counts
count_data_G2 <- data.frame(
  Category = c("Lost Regions", "Gained Regions"),
  Count = c(lost_regions_G2, gained_regions_G2)
)

# Reorder the factor levels of Category to have "Lost Regions" on the left
count_data_G2 <- count_data_G2 %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

plot_G2 <- ggplot(count_data_G2, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of", TF, "Regions in G2/M Phase"), x = "", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue"))

# Display the plot
print(plot_G2)

ggsave(filename = paste0(file_path, TF, "_Regions_G2.pdf"), plot = plot_G2)




# Combine count data for all phases
combined_count_data <- rbind(count_data_G1, count_data_S, count_data_G2)

# Add a new column for phase
combined_count_data$Phase <- c(rep("G0/G1 Phase", nrow(count_data_G1)),
                               rep("S Phase", nrow(count_data_S)),
                               rep("G2/M Phase", nrow(count_data_G2)))

# Reorder the factor levels of Phase
combined_count_data$Phase <- factor(combined_count_data$Phase, levels = c("G0/G1 Phase", "S Phase", "G2/M Phase"))

# Reorder the factor levels of Category within each phase
combined_count_data <- combined_count_data %>%
  mutate(Category = factor(Category, levels = c("Lost Regions", "Gained Regions")))

# Create the side-by-side bar chart with facets
plot <- ggplot(combined_count_data, aes(x = Phase, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Number of Differential", TF, "Binding"), x = "Cell Cycle Phase", y = "Number of Regions", fill = "Differential Binding") +
  theme_minimal() +
  scale_fill_manual(values = c("Lost Regions" = "darkgrey", "Gained Regions" = "lightblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))

plot

# Save the plot as a PDF
ggsave(filename = paste0(file_path, TF, "_Regions_Cell_Cycle_Phases.pdf"), plot = plot)
```
