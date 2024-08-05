# Data Preprocessing
## Alignment
```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=44gb
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

module load bowtie2/2.4.1
ref='/projects/ag-haensel/Pascal/genome_files/mm39_bowtie/ref/ref'

for f1 in *_R1_001.fastq.gz; do
f2=${f1%_R1_001.fastq.gz}_R2_001.fastq.gz
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 16 -I 10 -X 700 -x "$ref" -1 "$f1" -2 "$f2" -S "/scratch/phunold/TFCT/ESC_d2EpiLC_ONSMY/alignment/sam/${f1%%}_bowtie2.sam"
done
```
## Generate BAM files
```bash
for f1 in *_bowtie2.sam; do
sbatch -J StoB --mem 8GB --wrap "module load samtools/1.13 && samtools view -bS -F 0x04 "$f1" > /scratch/phunold/TFCT/ESC_d2EpiLC_ONSMY/alignment/bam/${f1%%}_bowtie2.mapped.bam"
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
    if [[ "$file" == *"ATAC"* ]]; then
        desired_reads=15000000
    elif [[  "$file" == *"MYC"* || "$file" == *"YAP1"* ]]; then
        desired_reads=5000000
    elif [[ "$file" == *"OCT4"* ||  ]]; then
        desired_reads=12000000
    elif [[ "$file" == *"SOX2"* ||  ]]; then
        desired_reads=19000000
    elif [[ "$file" == *"H3K27me3"* || "$file" == *"H3K4me3"* || "$file" == *"CTCF"* || "$file" == *"NANOG"* ]]; then
        desired_reads=10000000
    fi
    
    output_file="${file%.*}_${desired_reads}_down_bowtie2.mapped.bam"
    
    sbatch --mem 8G --cpus-per-task 8 --wrap "module load samtools/1.13 && total_reads=\$(samtools view -@ 8 -c -F 260 $file); scaling_factor=\$(bc <<< \"scale=4; $desired_reads / \$total_reads\"); samtools view -@ 8 -bs \"\$scaling_factor\" $file > $output_file"
done
```
## Sort BAM files
```bash
for f in *down_bowtie2.mapped.bam; do
 sbatch --mem 8G -J samSort --cpus-per-task 8 --wrap "module load samtools/1.13 && samtools sort -@ 8 $f > ${f%%.bam}.sort.bam"
done
```
## Call Peaks with MACS2
```bash
for f in *.sorted.bam; do
sbatch -J MACS2 --mem 32GB --wrap "module load use.own && module load pypack/macs2 && macs2 callpeak -t $f -f BAMPE -g mm --keep-dup all -n /scratch/phunold/ESC/peaks/$f --nomodel --extsize 55 -B --SPMR"
done

module load bedtools/2.29.2
for f in *_peaks.narrowPeak; do
  base_name=$(basename "$f" ".sorted.bam_peaks.narrowPeak")
  awk '{print $1"\t"$2"\t"$3}' "$f" | sortBed -i - | mergeBed -i - > "${base_name}_peaks.bed"
done
```
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
## Filter Consensus Peaks
```bash
module load bedtools/2.31.0

multiIntersectBed -i ESC-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G1_master_peaks.bed

multiIntersectBed -i ESC-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-S_master_peaks.bed

multiIntersectBed -i ESC-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-OCT4-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-OCT4-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-OCT4-G2_master_peaks.bed


multiIntersectBed -i ESC-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G1_master_peaks.bed

multiIntersectBed -i ESC-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-S_master_peaks.bed

multiIntersectBed -i ESC-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-SOX2-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-SOX2-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-SOX2-G2_master_peaks.bed


multiIntersectBed -i ESC-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G1_master_peaks.bed

multiIntersectBed -i ESC-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-S_master_peaks.bed

multiIntersectBed -i ESC-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-NANOG-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-NANOG-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-NANOG-G2_master_peaks.bed


multiIntersectBed -i ESC-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G1_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G1_master_peaks.bed

multiIntersectBed -i ESC-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-S_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-S_master_peaks.bed

multiIntersectBed -i ESC-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-YAP1-G2_master_peaks.bed
multiIntersectBed -i EpiLC-d2-YAP1-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-YAP1-G2_master_peaks.bed


multiIntersectBed -i ESC-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G1_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G1*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G1_master_peaks.bed

multiIntersectBed -i ESC-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-S_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-S*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-S_master_peaks.bed

multiIntersectBed -i ESC-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-MYC-G2_master_peaks.bed
multiIntersectBed -i EpiLC-2d-MYC-G2*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-MYC-G2_master_peaks.bed

multiIntersectBed -i ESC-ATAC*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ESC-ATAC_master_peaks.bed
multiIntersectBed -i EpiLC-ATAC*_peaks.bed | awk '{if($4>=2) print $0}' | sortBed -i - | mergeBed -i - > EpiLC-ATAC_master_peaks.bed


cat ESC-OCT4-G1_master_peaks.bed EpiLC-OCT4-G1_master_peaks.bed | sortBed -i - > OCT4-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G1_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-S_master_peaks.bed EpiLC-OCT4-S_master_peaks.bed | sortBed -i - > OCT4-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-S_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-S_all_peaks.over59nt.sorted.bed

cat ESC-OCT4-G2_master_peaks.bed EpiLC-OCT4-G2_master_peaks.bed | sortBed -i - > OCT4-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' OCT4-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > OCT4-G2_all_peaks.over59nt.sorted.bed


cat ESC-SOX2-G1_master_peaks.bed EpiLC-SOX2-G1_master_peaks.bed | sortBed -i - > SOX2-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G1_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-S_master_peaks.bed EpiLC-SOX2-S_master_peaks.bed | sortBed -i - > SOX2-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-S_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-S_all_peaks.over59nt.sorted.bed

cat ESC-SOX2-G2_master_peaks.bed EpiLC-SOX2-G2_master_peaks.bed | sortBed -i - > SOX2-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' SOX2-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > SOX2-G2_all_peaks.over59nt.sorted.bed


cat ESC-NANOG-G1_master_peaks.bed EpiLC-NANOG-G1_master_peaks.bed | sortBed -i - > NANOG-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G1_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-S_master_peaks.bed EpiLC-NANOG-S_master_peaks.bed | sortBed -i - > NANOG-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-S_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-S_all_peaks.over59nt.sorted.bed

cat ESC-NANOG-G2_master_peaks.bed EpiLC-NANOG-G2_master_peaks.bed | sortBed -i - > NANOG-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' NANOG-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > NANOG-G2_all_peaks.over59nt.sorted.bed


cat ESC-YAP1-G1_master_peaks.bed EpiLC-YAP1-G1_master_peaks.bed | sortBed -i - > YAP1-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G1_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-S_master_peaks.bed EpiLC-YAP1-S_master_peaks.bed | sortBed -i - > YAP1-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-S_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-S_all_peaks.over59nt.sorted.bed

cat ESC-YAP1-G2_master_peaks.bed EpiLC-YAP1-G2_master_peaks.bed | sortBed -i - > YAP1-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' YAP1-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > YAP1-G2_all_peaks.over59nt.sorted.bed


cat ESC-MYC-G1_master_peaks.bed EpiLC-MYC-G1_master_peaks.bed | sortBed -i - > MYC-G1_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G1_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G1_all_peaks.over59nt.sorted.bed

cat ESC-MYC-S_master_peaks.bed EpiLC-MYC-S_master_peaks.bed | sortBed -i - > MYC-S_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-S_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-S_all_peaks.over59nt.sorted.bed

cat ESC-MYC-G2_master_peaks.bed EpiLC-MYC-G2_master_peaks.bed | sortBed -i - > MYC-G2_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' MYC-G2_all_peaks.bed | sortBed -i - | mergeBed -i - > MYC-G2_all_peaks.over59nt.sorted.bed


cat ESC-ATAC_master_peaks.bed EpiLC-ATAC_master_peaks.bed | sortBed -i - > ATAC_all_peaks.bed
awk '{if (($3-$2) >= 60) print $0}' ATAC_all_peaks.bed | sortBed -i - | mergeBed -i - > ATAC_all_peaks.over59nt.sorted.bed
```
## Generate Count Matrices
```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

cell_lines=("ESC" "EpiLC-d2")
epitopes=("MYC" "SOX2" "NANOG" "OCT4" "YAP1")
phases=("G1" "G2" "S")
peak_dir="/scratch/phunold/ESC_EpiLC/bam/peaks/consensus/master_peaks"
bam_dir="/scratch/phunold/ESC_EpiLC/bam"
bedgraph_path="/scratch/phunold/ESC_EpiLC/bedgraph"

module load bedtools/2.29.2

for cell_line in "${cell_lines[@]}"; do
    for epitope in "${epitopes[@]}"; do
        for phase in "${phases[@]}"; do
            peak_file="${peak_dir}/${epitope}-${phase}_all_peaks.over59nt.sorted.bed"
            if [ -f "$peak_file" ]; then
                for f1 in "$bam_dir"/*"${cell_line}-${epitope}-${phase}"*.sorted.bam; do
                    bedtools coverage -a "$peak_file" -b "$f1" -counts > "$bedgraph_path/$(basename ${f1%%.sorted.bam}).bedgraph"
                    awk -v OFS='\t' '{print "chr"$1":"$2"-"$3, $4}' "$bedgraph_path/$(basename ${f1%%.sorted.bam}).bedgraph" > "$bedgraph_path/$(basename ${f1%%.sorted.bam})_counts.txt"
                done
            fi
        done
    done
done
```
## Generate BigWig Files for optical inspection
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
## plotCorrelation
```bash

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-CTCF-W300.merged.sorted.bam ESC-CTCF-ICS.merged.sorted.bam -p 8 -o CTCF_W300vsICS.npz "
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-H3K27me3-W300.merged.sorted.bam ESC-H3K27me3-ICS.merged.sorted.bam -p 8 -o H3K27me3_W300vsICS.npz "
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-H3K4me3-W300.merged.sorted.bam ESC-H3K4me3-ICS.merged.sorted.bam -p 8 -o H3K4me3_W300vsICS.npz "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in CTCF_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o CTCF_W300vsICS_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in H3K27me3_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o H3K27me3_W300vsICS_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in H3K4me3_W300vsICS.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o H3K4me3_W300vsICS_plotCorrelation.pdf "


sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-G1-1_S1_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-G1-2_S4_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-S-1_S2_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-S-2_S5_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-OCT4-G2-1_S3_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-OCT4-G2-2_S6_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_OCT4_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-G1-1_S7_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-G1-2_S10_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-S-1_S8_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-S-2_S11_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-NANOG-G2-1_S9_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-NANOG-G2-2_S12_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_NANOG_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-G1-1_S13_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-G1-2_S16_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-S-1_S14_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-S-2_S17_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-SOX2-G2-1_S15_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-SOX2-G2-2_S18_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_SOX2_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-G1-1_S1_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-G1-2_S4_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-S-1_S2_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-S-2_S5_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-MYC-G2-1_S3_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-MYC-G2-2_S6_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_MYC_G2.npz"

sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-G1-1_S25_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-G1-2_S28_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_G1.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-S-1_S26_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-S-2_S29_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_S.npz"
sbatch --mem 16G -J mBaSum --cpus-per-task 8 --wrap "module load use.own && module load pypack/deeptools && multiBamSummary bins -bs 10 -b ESC-YAP1-G2-1_S27_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam ESC-YAP1-G2-2_S30_L001_R1_001.fastq.gz_bowtie2.sam_bowtie2.mapped_down_bowtie2.mapped.bam.sorted.bam -p 8 -o ESC_YAP1_G2.npz"

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_OCT4_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_OCT4_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_NANOG_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_NANOG_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_SOX2_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_SOX2_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_MYC_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_MYC_G2_plotCorrelation.pdf "

sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_G1.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_G1_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_S.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_S_plotCorrelation.pdf "
sbatch --mem 64G -J pCorr --wrap "module load use.own && module load pypack/deeptools && plotCorrelation -in ESC_YAP1_G2.npz --corMethod pearson --skipZeros --whatToPlot scatterplot --log1p -o ESC_YAP1_G2_plotCorrelation.pdf "
```
## Inter Size Plot
```bash
for f in *ICS*merged.sorted.bam; do
 sbatch --mem 8G --time=2:00:00 -J size --wrap "module load openjdk/1.8.0_60 && java -Xmx7g -jar /home/phunold/privatemodules/pypack/picard.jar CollectInsertSizeMetrics \
      I=$f \
      O=./${f%%}_size_metrics.txt \
      H=./${f%%}_size_histogram.pdf"
done

for f in *W300*merged.sorted.bam; do
 sbatch --mem 8G --time=2:00:00 -J size --wrap "module load openjdk/1.8.0_60 && java -Xmx7g -jar /home/phunold/privatemodules/pypack/picard.jar CollectInsertSizeMetrics \
      I=$f \
      O=./${f%%}_size_metrics.txt \
      H=./${f%%}_size_histogram.pdf"
done
```
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
# MEME-ChIP
## Generate fasta Files
```bash
module load bedtools/2.29.2

for f in *.over59nt.sorted.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "module load bedtools/2.31.0 && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa -bed $f -fo ./${f%.conservedPeaks_MACS2_real.bed}.fasta"
done
```
# TOBIAS Analysis
```bash
mkdir ATACorrect_ATAC_EpiLC

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# load own modules
module load use.own

# load TOBIAS
module load pypack/tobias

# run TOBIAS
tobias ATACorrect --bam EpiLC-ATAC.merged.sorted.bam --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir ./ATACorrect_ATAC_EpiLC --cores 8


# deactivate conda env
conda deactivate

---
mkdir ATACorrect_ATAC_ESC

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# load own modules
module load use.own

# load TOBIAS
module load pypack/tobias

# run TOBIAS
tobias ATACorrect --bam ESC-ATAC.merged.sorted.bam --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir ./ATACorrect_ATAC_ESC --cores 8


# deactivate conda env
conda deactivate


#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module purge
module load miniconda/py38_4.9.2
module load use.own
module load pypack/tobias
for f in *_corrected.bw
do
tobias FootprintScores --signal $f --regions /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --output ${f%%_corrected.bw}_footprints.bw --cores 8
done

mkdir BINDetect_output

#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=phunold@uni-koeln.de
module purge
module load miniconda/py38_4.9.2
module load use.own
module load pypack/tobias
tobias BINDetect --motifs /projects/ag-haensel/CUT_Tag/Gisela/tobias/motifs/JASPAR_all_motifs/all_motifs.jaspar --signals ESC-ATAC.merged.sorted_footprints.bw EpiLC-ATAC.merged.sorted_footprints.bw --genome /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa --peaks /scratch/phunold/ESC/peaks/ATAC_all_peaks.over59nt.sorted.bed --outdir BINDetect_output --cond_names ESC EpiLC --cores 8
```














