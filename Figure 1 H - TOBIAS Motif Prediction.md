# Preprocessing in DynaTag_bulk_ESC_EpiLC.md
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

