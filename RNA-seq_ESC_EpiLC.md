# Data processing and Analyses of bulk RNA-seq of ESC and EpiLC
## Trimming
```bash
for f in *.fastq.gz; do
sbatch --mem 40GB -J trimming --cpus-per-task 15 --mail-type=FAIL,END --mail-user=phunold@uni-koeln.de --wrap "module load openjdk/11.0.2 && /projects/ag-haensel/tools/bbmap/bbduk.sh in=$f out=${f%%L000_R1_001.fastq.gz}clean.fastq.gz ref=/projects/ag-haensel/tools/bbmap/resources/polyA.fa.gz,/projects/ag-haensel/tools/bbmap/resources/truseq.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20"
done
```
