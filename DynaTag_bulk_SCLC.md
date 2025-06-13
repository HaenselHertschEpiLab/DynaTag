# Data Preprocessing
## Alignment
```bash
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module load bowtie2/2.4.1

ref='/projects/ag-haensel/Pascal/genome_files/hg38_bowtie/ref/ref'

for f1 in *_R1.fastq.gz; do
f2=${f1%_R1.fastq.gz}_R2.fastq.gz
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 -x "$ref" -1 "$f1" -2 "$f2" -S "/scratch/phunold/PDX_SCLC/fastq/alignment/sam/${f1%%}.sam"
done
```
## Generating BAM files
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 8GB --wrap "module load samtools/1.13 && samtools view -bS -F 0x04 "$f1" > /scratch/phunold/SCLC_PDX/alignment/bam/${f1%%}_bowtie2.mapped.bam"
done
```
## Downsampling
```bash
for f1 in *.bam; do
sbatch --mem 8G -J reads --mail-type=FAIL --mail-user=phunold@uni-koeln.de --wrap "module load samtools/1.13 && samtools view -c -F 260 $f1 >  ${f1%%.bam}.stat5"
done

cat *stat5 > reads.txt

#!/bin/bash
for file in *.bam; do
    prefix=$(echo "$file" | awk -F'[-_]' '{print $1}')
    condition=$(echo "$file" | awk -F'[-_]' '{print $2}')
    run=$(echo "$file" | awk -F'[-_]' '{print $3}')
    case "$prefix" in
        "ASCL1" | "FOXA1" | "MYC" | "NEUROD1" | "NFIB" | "NRF1" | "P53" | "P73" | "POU2F3" | "SP2" | "YAP1") desired_reads=10000000 ;;
    esac
    output_file="${prefix}_${condition}_${run}_${desired_reads}M.bam"
    sbatch --mem 8G --cpus-per-task 8 --wrap "module load samtools/1.13 && total_reads=\$(samtools view -@ 8 -c -F 260 $file); scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\"); samtools view -@ 8 -bs \"\$scaling_factor\" $file > $output_file"
done
```
## Sorting BAM files
```bash
for f in *M.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load samtools/1.13 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
```
## Peak Calling with MACS2
```bash
for f in *M.sort.bam; do
    sbatch -J MACS2 --mem 32GB --wrap "module load use.own && module load pypack/macs2 && macs2 callpeak -t $f -f BAMPE -g hs --keep-dup all -n /scratch/phunold/PDX_SCLC/fastq/alignment/peaks/$f --nomodel --extsize 55 -B --SPMR"
done

module load bedtools/2.29.2
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" "_10000000M.sort.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
```
## Calculate FRiP Scores
```bash
#!/bin/bash

bam_dir="/scratch/phunold/PDX_SCLC/fastq/alignment/bam/bam_down"
peak_dir="/scratch/phunold/PDX_SCLC/fastq/alignment/peaks"
output_file="FRiP_scores.txt"

module load samtools/1.13
module load bedtools/2.31.0

> "$output_file"

for bam_file in "$bam_dir"/*.bam; do
    sample_name=$(basename "$bam_file" .bam)
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

## Filter for consensus Peaks
```bash
module load bedtools/2.31.0

multiIntersectBed -i ASCL1_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ASCL1_CHEM_master_peaks.bed
multiIntersectBed -i ASCL1_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ASCL1_CTRL_master_peaks.bed

multiIntersectBed -i NEUROD1_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NEUROD1_CHEM_master_peaks.bed
multiIntersectBed -i NEUROD1_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NEUROD1_CTRL_master_peaks.bed

multiIntersectBed -i POU2F3_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > POU2F3_CHEM_master_peaks.bed
multiIntersectBed -i POU2F3_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > POU2F3_CTRL_master_peaks.bed

multiIntersectBed -i YAP1_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > YAP1_CHEM_master_peaks.bed
multiIntersectBed -i YAP1_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > YAP1_CTRL_master_peaks.bed

multiIntersectBed -i MYC_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > MYC_CHEM_master_peaks.bed
multiIntersectBed -i MYC_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > MYC_CTRL_master_peaks.bed

multiIntersectBed -i FOXA1_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > FOXA1_CHEM_master_peaks.bed
multiIntersectBed -i FOXA1_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > FOXA1_CTRL_master_peaks.bed

multiIntersectBed -i NFIB_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NFIB_CHEM_master_peaks.bed
multiIntersectBed -i NFIB_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NFIB_CTRL_master_peaks.bed

multiIntersectBed -i NRF1_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NRF1_CHEM_master_peaks.bed
multiIntersectBed -i NRF1_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > NRF1_CTRL_master_peaks.bed

multiIntersectBed -i P53_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > P53_CHEM_master_peaks.bed
multiIntersectBed -i P53_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > P53_CTRL_master_peaks.bed

multiIntersectBed -i P73_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > P73_CHEM_master_peaks.bed
multiIntersectBed -i P73_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > P73_CTRL_master_peaks.bed

multiIntersectBed -i SP2_CHEM*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > SP2_CHEM_master_peaks.bed
multiIntersectBed -i SP2_CTRL*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > SP2_CTRL_master_peaks.bed

cat ASCL1_CHEM_master_peaks.bed ASCL1_CTRL_master_peaks.bed | sortBed -i - > ASCL1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' ASCL1_all_peaks.bed | sortBed -i - | mergeBed -i - > ASCL1_all_peaks.over59nt.sorted.bed

cat NEUROD1_CHEM_master_peaks.bed NEUROD1_CTRL_master_peaks.bed | sortBed -i - > NEUROD1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NEUROD1_all_peaks.bed | sortBed -i - | mergeBed -i - > NEUROD1_all_peaks.over59nt.sorted.bed

cat POU2F3_CHEM_master_peaks.bed POU2F3_CTRL_master_peaks.bed | sortBed -i - > POU2F3_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' POU2F3_all_peaks.bed | sortBed -i - | mergeBed -i - > POU2F3_all_peaks.over59nt.sorted.bed

cat YAP1_CHEM_master_peaks.bed YAP1_CTRL_master_peaks.bed | sortBed -i - > YAP1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1_all_peaks.over59nt.sorted.bed

cat MYC_CHEM_master_peaks.bed MYC_CTRL_master_peaks.bed | sortBed -i - > MYC_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC_all_peaks.over59nt.sorted.bed

cat FOXA1_CHEM_master_peaks.bed FOXA1_CTRL_master_peaks.bed | sortBed -i - > FOXA1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' FOXA1_all_peaks.bed | sortBed -i - | mergeBed -i - > FOXA1_all_peaks.over59nt.sorted.bed

cat NFIB_CHEM_master_peaks.bed NFIB_CTRL_master_peaks.bed | sortBed -i - > NFIB_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NFIB_all_peaks.bed | sortBed -i - | mergeBed -i - > NFIB_all_peaks.over59nt.sorted.bed

cat NRF1_CHEM_master_peaks.bed NRF1_CTRL_master_peaks.bed | sortBed -i - > NRF1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NRF1_all_peaks.bed | sortBed -i - | mergeBed -i - > NRF1_all_peaks.over59nt.sorted.bed

cat P53_CHEM_master_peaks.bed P53_CTRL_master_peaks.bed | sortBed -i - > P53_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' P53_all_peaks.bed | sortBed -i - | mergeBed -i - > P53_all_peaks.over59nt.sorted.bed

cat P73_CHEM_master_peaks.bed P73_CTRL_master_peaks.bed | sortBed -i - > P73_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' P73_all_peaks.bed | sortBed -i - | mergeBed -i - > P73_all_peaks.over59nt.sorted.bed

cat SP2_CHEM_master_peaks.bed SP2_CTRL_master_peaks.bed | sortBed -i - > SP2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SP2_all_peaks.bed | sortBed -i - | mergeBed -i - > SP2_all_peaks.over59nt.sorted.bed
```
## Generate Count Matrices
```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=ASCL1
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=NEUROD1
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=POU2F3
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=YAP1
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=MYC
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=FOXA1
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=NFIB
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=NRF1
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=P53
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=P73
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
epitope=SP2
peak_path=/scratch/phunold/SCLC_PDX/alignment/peaks/${epitope}_all_peaks.over59nt.sorted.bed
bedgraph_path=/scratch/phunold/SCLC_PDX/alignment/bedgraph
bam_path=/scratch/phunold/SCLC_PDX/alignment/bam/down
module load bedtools/2.29.2
for f1 in $bam_path/${epitope}*_10000000M.sort.bam; do 
bedtools coverage -a $peak_path -b $f1 -counts > $bedgraph_path/$(basename ${f1%%_10000000M.sort.bam}).bedgraph
done

mkdir counts

for f1 in *.bedgraph; do
awk '{print $1 ":" $2 "-" $3, $4}' $f1 > ./counts/${f1%%.bedgraph}_counts.txt
done
```
# Proceed with DiffBind (?)

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


# PlotProfile for GSEA-Enriched Pathways
```bash
 Generate PlotProfiles for example TF/Pathway combinations
#!/bin/bash

bed_file="/projects/ag-haensel/Pascal/genome_files/hg38/genes_hg38_both.sort.bed"
input_file="/scratch/phunold/PDX_SCLC/fastq/alignment/GSEA_BED/TSS_bedgraph/POU2F3_e2f_targets_TSS.txt"

temp_file="temp.bed"

while read -r gene_name _; do
    grep -w "$gene_name" "$bed_file" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}'
done < "$input_file" > "$temp_file"

output_file="POU2F3_e2f_targets.bed"
cp "$temp_file" "$output_file"
rm "$temp_file"
done


module load bedtools/2.29.2
for f in *.bed; do
bedtools sort -i $f > ${f%%.bed}_sorted.bed
done


#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=15
module load use.own && module load pypack/deeptools
computeMatrix reference-point -S POU2F3.CTRL.PDXS02730.cpm.bw \
                              POU2F3.CHEM.PDXS02730.cpm.bw \
                              -R /scratch/phunold/PDX_SCLC/fastq/alignment/bam/bam_down/PP/POU2F3_e2f_targets_sorted.bed \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 250 \
                              --referencePoint TSS \
                              --missingDataAsZero \
                              --skipZeros -o /scratch/phunold/PDX_SCLC/fastq/alignment/bam/bam_down/PP/POU2F3_e2f_targets.mat.gz -p 15

sbatch -J pHM --mem 16GB --wrap "module load use.own && module load pypack/deeptools && plotProfile -m POU2F3_e2f_targets.mat.gz -out POU2F3_e2f_targets_profile.pdf --perGroup --colors blue red --plotTitle "POU2F3_e2f_targets" --plotType=fill"
```











