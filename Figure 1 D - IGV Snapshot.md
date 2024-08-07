# Preprocessing in DynaTag_bulk_ESC_EpiLC.md
# Generate BigWig Files for optical inspection
```bash

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-2d-MYC-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-2d-MYC-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-MYC-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-MYC-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-OCT4-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-OCT4-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-OCT4-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-OCT4-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-SOX2-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-SOX2-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-SOX2-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-SOX2-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-NANOG-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-NANOG-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-NANOG-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-NANOG-G2.merged.sorted.bam -"


sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-d2-YAP1-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o EpiLC-d2-YAP1-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-G1-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-G1.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-S-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-S.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-YAP1-G2-*down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-YAP1-G2.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - EpiLC-ATAC-*_down_bowtie2.mapped.sort.bam | samtools sort -o EpiLC-ATAC.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-ATAC*_down_bowtie2.mapped.sort.bam | samtools sort -o ESC-ATAC.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-CTCF-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-CTCF-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-CTCF-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-CTCF-W300.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K27me3-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K27me3-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K27me3-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K27me3-W300.merged.sorted.bam -"

sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K4me3-ICS*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K4me3-ICS.merged.sorted.bam -"
sbatch -J merge_sort --mem 8GB --wrap "module load samtools/1.13 && samtools merge - ESC-H3K4me3-W300*_down_bowtie2.mapped.bam.sorted.bam | samtools sort -o ESC-H3K4me3-W300.merged.sorted.bam -"

for f1 in *merged.sorted.bam; do
sbatch -J BAMindex --mem 8GB --wrap "module load samtools/1.13 && samtools index $f1"
done

for f1 in *merged.sorted.bam; do
sbatch --mem 16G -J BAM2BW --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && bamCoverage -p 8 -b $f1 -bs 10 --skipNAs --centerReads --normalizeUsing CPM -of bigwig -o ./${f1%%.sorted.bam}_cpm.bw "
done
```
