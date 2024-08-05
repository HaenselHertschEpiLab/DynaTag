# Data processing and Analyses of bulk RNA-seq of ESC and EpiLC
## Trimming
```bash
for f in *.fastq.gz; do
sbatch --mem 40GB -J trimming --cpus-per-task 15 --mail-type=FAIL,END --mail-user=phunold@uni-koeln.de --wrap "module load openjdk/11.0.2 && /projects/ag-haensel/tools/bbmap/bbduk.sh in=$f out=${f%%L000_R1_001.fastq.gz}clean.fastq.gz ref=/projects/ag-haensel/tools/bbmap/resources/polyA.fa.gz,/projects/ag-haensel/tools/bbmap/resources/truseq.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20"
done
```
## Alignment via STAR
```bash
#!/bin/bash
    if [ ! -d "STAR_aligned_bam" ]; then
      mkdir STAR_aligned_bam
    fi
    for f in ./*clean.fastq.gz
    do
      module load star/2.7.8a
      STAR --runThreadN 8 \
      --genomeDir /projects/ag-haensel/Michaela/genome_files/mm10_STAR \
      --readFilesIn $f \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.6 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outSAMattributes NH HI NM MD \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix /scratch/phunold/RPM/STAR_aligned_bam_mm10/$f \
      --readFilesCommand gunzip -c
    done
```
