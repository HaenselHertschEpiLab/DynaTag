# Preprocessing in DynaTag_bulk_SCLC.md
# Bubble Plot for GSEA-Enriched PAthways
## Plot GSEA and generate gene lists as txt
```bash
# Plot GSEA and generate gene lists as txt
CTRL <- read_tsv("/Users/pascalhunold/gsea_home/output/apr22/PDX_CTRL_CHEM_GSEA.GseaPreranked.1713793795039/gsea_report_for_na_neg_1713793795039.tsv")
CHEM <- read_tsv("/Users/pascalhunold/gsea_home/output/apr22/PDX_CTRL_CHEM_GSEA.GseaPreranked.1713793795039/gsea_report_for_na_pos_1713793795039.tsv")

CTRL$group <- "CTRL"
CHEM$group <- "CHEM"

combined_data <- rbind(CTRL, CHEM)
combined_data <- combined_data %>% filter(`NOM p-val` <= 0.05)
combined_data$NAME <- with(combined_data, reorder(NAME, NES))
combined_data$size <- 1 / ifelse(combined_data$`FDR q-val` == 0, min(combined_data$`FDR q-val`[combined_data$`FDR q-val` > 0]), combined_data$`FDR q-val`)
max_NES <- max(abs(combined_data$NES))

p1 <- ggplot(combined_data, aes(x = reorder(NAME, NES), y = NES, fill = group)) +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("CTRL" = "darkgrey", "CHEM" = "lightblue")) +
  labs(x = "", y = "Normalised Enrichment Score", title = "GSEA: CTRL vs CHEM") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = "black")) +
  scale_y_continuous(limits = c(-max_NES, max_NES), breaks = seq(floor(-max_NES), ceiling(max_NES), by = 1))

combined_data
ggsave("PDX_SCLC_GSEA_barchart.pdf", plot = p1)


###
read_pathway_files <- function(combined_data, directory_path) {
  for (pathway_name in combined_data$NAME) {

    file_name <- paste0(pathway_name, ".tsv")
    file_path <- file.path(directory_path, file_name)
    
    pathway_data <- read_tsv(file_path)
    
    pathway_data$Pathway <- pathway_name
    
    assign(paste0(tolower(gsub(" ", "_", pathway_name)), "_data"), pathway_data, envir = .GlobalEnv)
  }
}

directory_path <- "/Users/pascalhunold/gsea_home/output/apr22/PDX_CTRL_CHEM_GSEA.GseaPreranked.1713793795039/"

read_pathway_files(combined_data, directory_path)

###
generate_gene_lists <- function(data) {
  enriched_genes <- data$SYMBOL[data$`CORE ENRICHMENT` == "Yes"]
  
  return(enriched_genes)
}

gene_lists <- list()

for (pathway_data_name in ls(pattern = "^hallmark_.*_data$")) {
  pathway_data <- get(pathway_data_name)
  
  gene_list <- generate_gene_lists(pathway_data)
  
  gene_lists[[pathway_data_name]] <- gene_list
}

# Now gene_lists contains lists of enriched gene symbols for each pathway
# separate lists

read_pathway_files <- function(combined_data, directory_path) {
  for (pathway_name in combined_data$NAME) {
    file_name <- paste0(pathway_name, ".tsv")
    file_path <- file.path(directory_path, file_name)
    
    pathway_data <- read_tsv(file_path)
    
    enriched_genes <- pathway_data$SYMBOL[pathway_data$`CORE ENRICHMENT` == "Yes"]
    
    pathway_data$Pathway <- pathway_name
    
    assign(paste0(tolower(gsub(" ", "_", pathway_name)), "_filtered"), enriched_genes, envir = .GlobalEnv)
  }
}

directory_path <- "/Users/pascalhunold/gsea_home/output/apr22/PDX_CTRL_CHEM_GSEA.GseaPreranked.1713793795039/"

read_pathway_files(combined_data, directory_path)

### BED file
gtf_file <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/Homo_sapiens.GRCh38.109.gtf.gz"

gtf <- import.gff(gtf_file)

head(gtf)

output_directory <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/GSEA_BED/"

for (pathway_name in names(pathway_lists)) {
  pathway_genes <- pathway_lists[[pathway_name]]
  
  pathway_genes_gtf <- subset(gtf, type == "gene" & gene_name %in% pathway_genes)
  
  bed_data <- data.frame(
    chrom = seqnames(pathway_genes_gtf),
    start = start(pathway_genes_gtf),
    end = end(pathway_genes_gtf),
    name = mcols(pathway_genes_gtf)$gene_name
  )
  
  bed_file_path <- paste0(output_directory, pathway_name, ".bed")
  write.table(bed_data, bed_file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

### TXT file with gene symbols
output_directory <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/GSEA_BED/"

for (pathway_name in names(pathway_lists)) {
  pathway_genes <- pathway_lists[[pathway_name]]
  
  txt_file_path <- paste0(output_directory, pathway_name, ".txt")
  
  writeLines(as.character(pathway_genes), con = txt_file_path)
}
```
## Generate TSS-1000+250 txt file per GSEA pathway
```bash
#!/bin/bash

bed_file="/projects/ag-haensel/Pascal/genome_files/hg38/genes_hg38_both.sort.bed"

for filtered_file in *_filtered.txt; do
    temp_file="${filtered_file%_filtered.txt}_temp.txt"

    while read -r gene; do
        while read -r gene_name; do
            strand=$(grep -w "$gene_name" "$bed_file" | awk '{print $5}')
            if [ "$strand" == "+" ]; then
                grep -w "$gene_name" "$bed_file" | awk '{print $1, $2, $4}'
            elif [ "$strand" == "-" ]; then
                grep -w "$gene_name" "$bed_file" | awk '{print $1, $3, $4}'
            fi
	done < <(grep -w "$gene" "$filtered_file")
    done < "$filtered_file" > "$temp_file"

    output_file="${filtered_file%_filtered.txt}_TSS.txt"

    awk '{print $1, $2, $3, $2 - 1000, $2 + 250, $4}' "$temp_file" > "$output_file"

    rm "$temp_file"
done
```
## Generate sorfed BED files from the TSS.txt
```bash
#!/bin/bash

module load bedtools/2.29.2

for tss_file in *_TSS.txt; do
    output_bed="${tss_file%_TSS.txt}_sorted_TSS.bed"

    awk 'BEGIN{OFS="\t"} {print $1, $4, $5, $3}' "$tss_file" > "$output_bed"

    bedtools sort -i "$output_bed" > "${output_bed%.bed}_sorted.bed"
done
```
## Generate bedgraphs from the TSS.bed file per pathway 
```bash
#!/bin/bash

bam_path="/scratch/phunold/PDX_SCLC/fastq/alignment/bam"

bed_dir="/scratch/phunold/PDX_SCLC/fastq/alignment/GSEA_BED"

bedgraph_dir="/scratch/phunold/PDX_SCLC/fastq/alignment/GSEA_BED/TSS_bedgraph"

for bed_file in "${bed_dir}"/*_TSS_sorted_no_chr.bed; do
    bed_filename=$(basename "$bed_file" .bed)
    
    for bam_file in "${bam_path}"/*_merged.sort.bam; do
        bam_filename=$(basename "$bam_file" .bam)
        
        output_bedgraph="${bedgraph_dir}/${bam_filename}_${bed_filename}.bedgraph"
        
        sbatch -J StoB --mem 32GB --wrap "module load bedtools/2.29.2 && bedtools coverage -a \"$bed_file\" -b \"$bam_file\" -counts > \"$output_bedgraph\""
    done
done
```
## calculate ratio per gene per TF per pathway
```bash
#!/bin/bash

epitopes=("ASCL1" "NEUROD1" "POU2F3" "YAP1" "MYC" "FOXA1" "NFIB" "NRF1" "P53" "P73" "SP2")

hallmarks=("androgen_response" "angiogenesis" "apical_surface" "bile_acid_metabolism" "coagulation" "e2f_targets" "epithelial_mesenchymal_transition" "estrogen_response_early" "fatty_acid_metabolism" "glycolysis" "heme_metabolism" "hypoxia" "interferon_alpha_response" "mtorc1_signaling" "oxidative_phosphorylation" "tnfa_signaling_via_nfkb" "uv_response_dn")

for epitope in "${epitopes[@]}"; do
    for hallmark in "${hallmarks[@]}"; do
        chem_bedgraph="${epitope}_CHEM_merged.sort_hallmark_${hallmark}_sorted_TSS_sorted_no_chr.bedgraph"
        ctrl_bedgraph="${epitope}_CTRL_merged.sort_hallmark_${hallmark}_sorted_TSS_sorted_no_chr.bedgraph"

        output_file="${epitope}_${hallmark}_TSS.txt"

paste <(awk '{print $4, $5}' "$ctrl_bedgraph") \
      <(awk '{print $4, $5}' "$chem_bedgraph") \
| awk -v FS="[[:space:]]+" '{ if ($2 != 0 && $4 != 0) print $1, $4/$2; else print $1, "NA" }' \
| awk '{sum+=$2; count++} END {print "MEAN", sum/count} NR>1' > "$output_file"

    done
done

# mean_ratio.txt file
#!/bin/bash

epitopes=("ASCL1" "NEUROD1" "POU2F3" "YAP1" "MYC" "FOXA1" "NFIB" "NRF1" "P53" "P73" "SP2")

hallmarks=("androgen_response" "angiogenesis" "apical_surface" "bile_acid_metabolism" "coagulation" "e2f_targets" "epithelial_mesenchymal_transition" "estrogen_response_early" "fatty_acid_metabolism" "glycolysis" "heme_metabolism" "hypoxia" "interferon_alpha_response" "mtorc1_signaling" "oxidative_phosphorylation" "tnfa_signaling_via_nfkb" "uv_response_dn")

declare -A mean_ratios

for epitope in "${epitopes[@]}"; do
    for hallmark in "${hallmarks[@]}"; do
        input_file="${epitope}_${hallmark}_TSS.txt"

        mean_ratio=$(awk '/MEAN/ {print $2}' "$input_file")
        mean_ratios["$epitope,$hallmark"]=$mean_ratio
    done
done

output_file="mean_ratio.txt"
echo -e "Epitope\t${hallmarks[@]}" > "$output_file"  

for epitope in "${epitopes[@]}"; do
    row="$epitope"
    for hallmark in "${hallmarks[@]}"; do
        mean_ratio="${mean_ratios[$epitope,$hallmark]}"
        row+="\t${mean_ratio:-NA}" 
    done
    echo -e "$row" >> "$output_file"
done

echo "mean_ratio.txt generated successfully."
```
## t-test
```bash
#!/bin/bash

module load R/4.3.1_system

perform_ttest() {
    chem_bedgraph="$1"
    ctrl_bedgraph="$2"
    result_file="$3"

    genes=$(awk '$5 > 0' "$chem_bedgraph" | cut -f4)
    genes=$(comm -12 <(awk '$5 > 0' "$ctrl_bedgraph" | cut -f4 | sort) <(echo "$genes" | sort))

    awk -v genes="$genes" 'BEGIN { split(genes, arr, " "); for (i in arr) gene[arr[i]]=1 } $4 in gene { print }' "$chem_bedgraph" > chem_filtered.bedgraph
    awk -v genes="$genes" 'BEGIN { split(genes, arr, " "); for (i in arr) gene[arr[i]]=1 } $4 in gene { print }' "$ctrl_bedgraph" > ctrl_filtered.bedgraph

    Rscript - <<R_SCRIPT
        chem_data <- read.table("chem_filtered.bedgraph", header = FALSE)
        ctrl_data <- read.table("ctrl_filtered.bedgraph", header = FALSE)

        t_test_result <- t.test(chem_data\$V5, ctrl_data\$V5, paired = TRUE)

        write.table(t_test_result\$p.value, file = "$result_file", row.names = FALSE, col.names = FALSE, quote = FALSE)
R_SCRIPT

    rm chem_filtered.bedgraph ctrl_filtered.bedgraph
}

epitopes=("ASCL1" "NEUROD1" "POU2F3" "YAP1" "MYC" "FOXA1" "NFIB" "NRF1" "P53" "P73" "SP2")

hallmarks=("androgen_response" "angiogenesis" "apical_surface" "bile_acid_metabolism" "coagulation" "e2f_targets" "epithelial_mesenchymal_transition" "estrogen_response_early" "fatty_acid_metabolism" "glycolysis" "heme_metabolism" "hypoxia" "interferon_alpha_response" "mtorc1_signaling" "oxidative_phosphorylation" "tnfa_signaling_via_nfkb" "uv_response_dn")

for epitope in "${epitopes[@]}"; do
    for hallmark in "${hallmarks[@]}"; do
        chem_bedgraph="${epitope}_CHEM_merged.sort_hallmark_${hallmark}_sorted_TSS_sorted_no_chr.bedgraph"
        ctrl_bedgraph="${epitope}_CTRL_merged.sort_hallmark_${hallmark}_sorted_TSS_sorted_no_chr.bedgraph"

        pvalue_file="${epitope}_${hallmark}_ttest_pvalue.txt"

        perform_ttest "$chem_bedgraph" "$ctrl_bedgraph" "$pvalue_file"
    done
done
```
## Bubble Plot in R
```R
data <- read.table("/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration/DT_RNA/mean_ratio.txt", header=TRUE, row.names=1, sep="")
pvalue_file <- "/Users/pascalhunold/Desktop/PhD_Documentation/SCLC/snRNAseq/SCLC_PDX/Integration/DT_RNA/p_value_matrix.txt"

pvalues <- read.table(pvalue_file, header=TRUE, row.names=1, sep="")

new_data_list <- list()

# Iterate through each row (epitope) and column (hallmark) to populate the data
for (i in 1:nrow(data)) {
  epitope <- rownames(data)[i]
  for (j in 1:ncol(data)) {
    hallmark <- colnames(data)[j]
    ratio <- log2(data[i, j]) 
    pvalue <- -log10(pvalues[i, j]) 
    new_data_list[[length(new_data_list) + 1]] <- c(epitope, hallmark, ratio, pvalue)
  }
}

new_data_matrix <- do.call(rbind, new_data_list)

new_data <- as.data.frame(new_data_matrix)
colnames(new_data) <- c("epitope", "hallmark", "ratio", "pvalue")

new_data$pvalue <- as.numeric(new_data$pvalue)
new_data$ratio <- as.numeric(new_data$ratio)
print(new_data)

desired_order <- c("hypoxia", "epithelial_mesenchymal_transition", "apical_surface", "tnfa_signaling_via_nfkb", "mtorc1_signaling", "angiogenesis", "bile_acid_metabolism", "coagulation", "androgen_response", "uv_response_dn", "heme_metabolism", "estrogen_response_early", "fatty_acid_metabolism", "glycolysis", "e2f_targets", "interferon_alpha_response", "oxidative_phosphorylation")
new_data$hallmark <- factor(new_data$hallmark, levels = desired_order)

max_abs_ratio <- max(abs(new_data$ratio))
max_abs_pvalue <- max(abs(new_data$pvalue))

significant_data <- new_data %>%
  filter(pvalue >= 1.3)

bubble_plot <- ggplot(significant_data, aes(x = reorder(hallmark, -pvalue), y = reorder(epitope, ratio), 
                                            fill = ratio, size = pvalue)) +
  geom_point(shape = 21, color = "black") +  
  scale_fill_gradient2(low = "#0a12a8", mid = "white", high = "#bd0202", 
                       midpoint = 0, name = "Ratio [log2]", limits = c(-max_abs_ratio, max_abs_ratio)) + 
  labs(x = "", y = "", size = "-log10(p-value)", fill = "log2(Ratio)") +
  scale_size_continuous(name = "p-value [-log10]", range = c(1.3, 2*max_abs_pvalue)) +  # Size Scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )

bubble_plot

ggsave("/Users/pascalhunold/Desktop/PhD_Documentation/DynaTag/Sequencing/SCLC_PDX/GSEA_Plots/PDX_SCLC_GSEA_DynaTag_Integration.pdf", plot = bubble_plot, width = 14.72, height = 10.62)
```
