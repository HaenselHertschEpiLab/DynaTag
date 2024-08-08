# Proprocessing in DynaTag_bulk_ESC_EpiLC.md and DynaTag_sc_ESC.md
# Enrichment of ChIP-Atlas target genes
```R
# Read the TSV file into a data frame
df <- read.table("/Users/pascalhunold/Downloads/Nanog.1.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Extract the first column containing gene symbols
gene_symbols <- df$Target_genes

# Print or use gene_symbols as needed
print(gene_symbols)


# Path to your BED file
#mm39_bed <- "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/genomes/genes_mm39_both.sort.bed"
# Read BED file
#bed_data <- import.bed(mm39_bed)
# Convert to data frame for easier manipulation
#bed_df <- as.data.frame(bed_data)

# Filter BED data to keep only rows where gene symbol matches those in DESeq2 output
filtered_bed <- bed_df[bed_df$name %in% gene_symbols,]
# Select only the required columns in the order: chromosome, start, end, name, strand
output_bed <- filtered_bed[, c("seqnames", "start", "end", "name", "strand")]

# Write to a new BED file
write.table(output_bed, 
            file = "/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/ChIP_Atlas/NANOG_ChIP_Atlas_mm10_Targets.bed", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


TF <- "NANOG"

# Read the data from the BED file into a data frame
peaks_df <- read.table(paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/peaks/", TF, "-G1_all_peaks.over59nt.sorted.bed"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
peaks_df$V1 <- paste0("chr", peaks_df$V1)
peaks <- GRanges(seqnames=peaks_df$V1, ranges=IRanges(start=peaks_df$V2, end=peaks_df$V3), strand=peaks_df$V4)

target_genes <- fread(paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/ChIP_Atlas/", TF, "_ChIP_Atlas_mm10_Targets.bed"), header = FALSE)
target_genes_gr <- GRanges(seqnames = target_genes$V1, 
                           ranges = IRanges(start = target_genes$V2, end = target_genes$V3))



# Read the known genes data
known_genes <- fread("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/genomes/genes_mm39_both.sort.bed", header = FALSE)

# Convert to GenomicRanges object
known_genes_gr <- GRanges(seqnames = known_genes$V1, 
                          ranges = IRanges(start = known_genes$V2, end = known_genes$V3))

# Extract gene symbols from known_genes
known_gene_symbols <- known_genes$V4

# Extract gene symbols from target_genes
target_gene_symbols <- target_genes$V4

# Find gene symbols in known_genes that are not in target_genes
non_target_gene_symbols <- setdiff(known_gene_symbols, target_gene_symbols)

# Randomly sample gene symbols from non_target_gene_symbols
# Determine the number of rows in target_genes
num_rows_target_genes <- nrow(target_genes)

# Randomly sample gene symbols from non_target_gene_symbols
random_gene_symbols <- sample(non_target_gene_symbols, num_rows_target_genes, replace = FALSE)
length(random_gene_symbols)
# Filter known_genes based on random_gene_symbols
random_genes <- known_genes[V4 %in% random_gene_symbols]

# Convert filtered data to GRanges
random_genes_gr <- GRanges(
  seqnames = random_genes$V1,
  ranges = IRanges(start = random_genes$V2, end = random_genes$V3),
  strand = random_genes$V6,
  names = random_genes$V4
)



overlap_target <- countOverlaps(peaks, target_genes_gr)
overlap_random <- countOverlaps(peaks, random_genes_gr)
overlap_target_random <- countOverlaps(target_genes_gr, random_genes_gr)
length(peaks)
length(target_genes_gr)
length(random_genes_gr)
cat("Overlap with target genes:", sum(overlap_target > 0), "\n")
cat("Overlap with randomly selected genes:", sum(overlap_random > 0), "\n")
cat("Overlap between target and randomly selected genes:", sum(overlap_target_random > 0), "\n")
#library(ggplot2)

# Create a data frame for plotting
plot_data <- data.frame(
  Category = c("ChIP Atlas Target Genes", "Random Genes"),
  Count = c(sum(overlap_target > 0), sum(overlap_random > 0))
)

# Plot
gg <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), colour = "black") +
  labs(title = paste("Overlap between Peaks and Genes -", TF),
       x = "",
       y = "Peaks in Genes [#]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = c("Random Genes" = "#e3e3e3", "ChIP Atlas Target Genes" = "#666666")) +
  guides(fill = FALSE)
gg
# Save the plot as a PDF
ggsave(filename = paste0("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/ESC_EpiLC/figures/", TF, "_ChIPAtlas_enrichment.pdf"), plot = gg)
```
