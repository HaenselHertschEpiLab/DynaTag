# Data Preprocessing
## Demultiplexing via SpitPipe
```bash
# 1 - Demux Sub Libraries and align against hg38_mm39 (although called mm10 in the path)

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-1-cat_S1_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-1-cat_S1_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S1 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-2-cat_S2_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-2-cat_S2_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S2 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-3-cat_S3_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-3-cat_S3_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S3 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-4-cat_S4_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-4-cat_S4_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S4 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-5-cat_S5_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-5-cat_S5_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S5 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-6-cat_S6_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-6-cat_S6_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S6 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-7-cat_S7_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-7-cat_S7_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S7 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
--m all \
--chemistry V2 \
--fq1 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-8-cat_S8_L001_R1_001.fastq.gz \
--fq2 /projects/ag-haensel/raw_data/RNA_Seq/PARSE/2023-05-14_WT_full_scale_RPchmeo_CDXchemo_RPM/cat_fastq/SubLib-8-cat_S8_L001_R2_001.fastq.gz \
--output_dir /projects/ag-haensel/Pascal/SCLC/PDX/S8 \
--sample SCLC_CTRL_1 C9-C12 \
--sample SCLC_CTRL_2 D1-D4 \
--sample SCLC_CHEM_1 D5-D8 \
--sample SCLC_CHEM_2 D9-D12 \
--genome_dir /projects/ag-haensel/Pascal/genome_files/PARSE/genomes/hg38_mm10

#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --output=filename_%j.log   # Standard output and error log
#SBATCH --mem=64GB
#SBATCH --time=24:00:00

# clean module list
module purge

# load conda module
module load miniconda/py38_4.9.2

# activate new ParseBiosciences env
conda activate /opt/rrzk/software/conda-envs/ParseBioscience-1.0.6p/

# run split-pipe
# Achtung: hier besser Pfad zu sample "xcond_1" und "A1-A10" angeben
split-pipe \
    --mode comb \
    --sublibraries /projects/ag-haensel/Pascal/SCLC/PDX/S1 /projects/ag-haensel/Pascal/SCLC/PDX/S2 /projects/ag-haensel/Pascal/SCLC/PDX/S3 /projects/ag-haensel/Pascal/SCLC/PDX/S4 /projects/ag-haensel/Pascal/SCLC/PDX/S5 /projects/ag-haensel/Pascal/SCLC/PDX/S6 /projects/ag-haensel/Pascal/SCLC/PDX/S7 /projects/ag-haensel/Pascal/SCLC/PDX/S8 \
    --output_dir /projects/ag-haensel/Pascal/SCLC/PDX

# deactivate conda env
conda deactivate
```
