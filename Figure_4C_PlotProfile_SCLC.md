# Preprocessing in DynaTag_bulk_SCLC.md
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
