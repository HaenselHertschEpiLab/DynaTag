# Preprocessing in DynaTag_bulk_ESC_EpiLC.md
# Generate fasta Files
```bash
module load bedtools/2.29.2

for f in *.over59nt.sorted.bed; do
sbatch --mem 8G --time=2:00:00 -J FA --wrap "module load bedtools/2.31.0 && bedtools getfasta -fi /projects/ag-haensel/Pascal/genome_files/mm39_bowtie/Mus_musculus.GRCm39.dna.primary_assembly.fa -bed $f -fo ./${f%.conservedPeaks_MACS2_real.bed}.fasta"
done
```
