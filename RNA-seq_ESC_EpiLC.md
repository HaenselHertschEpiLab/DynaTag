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
      --outFileNamePrefix /scratch/phunold/$f \
      --readFilesCommand gunzip -c
    done
```
## Downsampling
```bash
for file in *.bam; do
sbatch --mem 8G -J stat5 --mail-type=FAIL --mail-user=phunold@uni-koeln.de --wrap "module load samtools/1.13 && samtools view -c -F 260 $file >  ${file%%.bam}.stat5"
done

if [ ! -d " downsampled
" ]; then
  mkdir downsampled
fi

for file in *.stat5
do
  stat5_file=`cat $file`
  echo $stat5_file
  if [ "$stat5_file" -gt 8000000 ]
  then
      echo "bigger than 7M" $file
      factor_down=$(awk -v m=$stat5_file 'BEGIN { print 7000000/m }')
      echo $factor_down
      echo "***"
      sbatch --mem 8G --cpus-per-task 8 --wrap "module load samtools/1.13 && samtools view -@ 8 -s $factor_down -b ${file%%.stat5}.bam > ./downsampled/${file%%.stat5}.8M.bam"
  else
    echo " =================> LESS"
    sbatch --mem 8G --wrap "cp ${file%%.stat5}.bam ${file%%.stat5}.less7M.bam"
  fi
  echo "================="
done
```
## Count Reads
```bash
if [ ! -d " counts
" ]; then
  mkdir counts
fi

for bam in *.bam
do
sbatch --mem 16g --mail-type=FAIL,END --mail-user=phunold@uni-koeln.de --wrap "module load use.own && module load pypack/htseq && htseq-count -f bam -m union -s no -t exon -i gene_id $bam /projects/ag-haensel/Michaela/genome_files/mm10_STAR/Mus_musculus.GRCm38.102_edited.gtf > ./counts/${bam%%.merged.bam}.counts.txt"
done
```




