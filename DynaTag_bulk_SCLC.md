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

for f1 in *_R1_001.trimmed.fastq.gz; do
f2=${f1%_R1_001.trimmed.fastq.gz}_R2_001.trimmed.fastq.gz
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












