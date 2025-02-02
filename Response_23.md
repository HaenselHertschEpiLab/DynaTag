## Code needs to be adjusted to SCLC peaks

## Response 23 HOMER2 known motif analyses SCLC TFs using matched DynaTag peaks and ChIP-Atlas derived peaks.
```bash
cat normalize_bed_files.sh
#!/bin/bash

# Parameters
input_files="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/peaks/ChIP_Atlas_target_peaks/*.hg38.merged.bed /scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/peaks/peaks_master/*.bed"  # Process all BED files in the directory
output_dir="normalized_peaks"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each BED file
for input_bed in $input_files; do
    # Extract the base name of the file (without path)
    base_name=$(basename "$input_bed")

    # Calculate average peak size
    avg_peak_size=$(awk '{sum += ($3 - $2); count++} END {print int(sum / count)}' "$input_bed")
    echo "Calculated average peak size for $input_bed: ${avg_peak_size} bp"

    # Use calculated average peak size for normalization
    normalized_bed="${output_dir}/normalized_${base_name}"
    awk -v size="$avg_peak_size" 'BEGIN {OFS="\t"} {
        mid = int(($2 + $3) / 2);
        start = mid - int(size / 2);
        end = mid + int(size / 2);
        if (start < 0) start = 0;
        print $1, start, end, ".", "1000";
    }' "$input_bed" > "$normalized_bed"

    echo "Normalized BED file created: $normalized_bed"
done
```
```bash
cat generate_background.sh
#!/bin/bash

# Parameters
normalized_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/peaks/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/Motifs_homer/bg_results"
genome_fa="/projects/ag-haensel/Pascal/genome_files/hg38_bowtie/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
email="rhaensel@uni-koeln.de"
conda_env="/projects/ag-haensel/tools/.conda/envs/abc-model-env"
kmer_value=2  # Fixed kmer value
time_limit="4:00:00"
memory_limit="32G"
cpus=8

# Create the background directory if it doesn't exist
mkdir -p "$background_dir"

# Loop through each normalized BED file
for normalized_bed in "$normalized_dir"/*.bed; do
    # Extract base name of the normalized file (without directory and extension)
    base_name=$(basename "$normalized_bed" .bed)
    bg_output_prefix="${background_dir}/bg_kmer${kmer_value}_${base_name}"

    # Prepare the SLURM script for background generation
    slurm_script="#!/bin/bash
#SBATCH --job-name=bg_kmer${kmer_value}_${base_name}
#SBATCH --output=${bg_output_prefix}_slurm.out
#SBATCH --error=${bg_output_prefix}_slurm.err
#SBATCH --time=${time_limit}
#SBATCH --mem=${memory_limit}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${email}

# Load and activate the conda environment
conda activate ${conda_env}

# Generate the background file
homer2 background \\
    -p ${normalized_bed} \\
    -g ${genome_fa} \\
    -pkmer ${kmer_value} \\
    -N 100000 \\
    -NN 100000000 \\
    -allowTargetOverlaps \\
    -allowBgOverlaps \\
    -o ${bg_output_prefix}

# Deactivate conda environment
conda deactivate
"

    # Save the SLURM script to a file
    slurm_script_file="submit_bg_kmer${kmer_value}_${base_name}.sbatch"
    echo "$slurm_script" > "$slurm_script_file"

    # Submit the SLURM job
    sbatch "$slurm_script_file"
done

echo "Background generation jobs have been submitted for all normalized peak files."
```
```bash
cat run_motif_analysis.sh
#!/bin/bash

# Parameters
normalized_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/peaks/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/Motifs_homer/bg_results"
output_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/Motifs_homer/motifs_results"
genome_fa="/projects/ag-haensel/Pascal/genome_files/hg38_bowtie/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
email="rhaensel@uni-koeln.de"
conda_env="/projects/ag-haensel/tools/.conda/envs/abc-model-env"
time_limit="4:00:00"
memory_limit="32G"
cpus=8
kmer_value=2  # Fixed kmer value

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each normalized BED file
for normalized_bed in "$normalized_dir"/*.bed; do
    # Extract base name of the normalized file (without directory and extension)
    base_name=$(basename "$normalized_bed" .bed)
    bg_file="${background_dir}/bg_kmer${kmer_value}_${base_name}.bg.positions.bed"
    motifs_output="${output_dir}/motifs_kmer${kmer_value}_${base_name}"

    # Check if the corresponding background file exists
    if [[ ! -f "$bg_file" ]]; then
        echo "Error: Background file $bg_file not found for $normalized_bed. Skipping..."
        continue
    fi

    # Prepare the SLURM script for motif analysis
    slurm_script="#!/bin/bash
#SBATCH --job-name=motifs_kmer${kmer_value}_${base_name}
#SBATCH --output=${motifs_output}_slurm.out
#SBATCH --error=${motifs_output}_slurm.err
#SBATCH --time=${time_limit}
#SBATCH --mem=${memory_limit}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${email}

# Load and activate the conda environment
conda activate ${conda_env}

# Run motif analysis (known motif enrichment only)
findMotifsGenome.pl \\
    ${normalized_bed} \\
    ${genome_fa} \\
    ${motifs_output} \\
    -size given \\
    -len 8,10,12 \\
    -S 25 \\
    -mis 2 \\
    -mset vertebrates \\
    -gc \\
    -bg ${bg_file} \\
    -nomotif \\
    -redundant 2.0 \\
    -maxN 0.7 \\
    -nlen 3 \\
    -nmax 160 \\
    -e 0.0 \\
    -minlp -10.0

# Deactivate conda environment
conda deactivate
"

    # Save the SLURM script to a file
    slurm_script_file="submit_motifs_kmer${kmer_value}_${base_name}.sbatch"
    echo "$slurm_script" > "$slurm_script_file"

    # Submit the SLURM job
    sbatch "$slurm_script_file"
done

echo "Motif analysis jobs have been submitted for all normalized peak files with corresponding background files."
```
```bash
nano run_denovo_motif_analysis.sh
#!/bin/bash

# Parameters
input_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/peaks/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/Motifs_homer/bg_results"
output_dir="/scratch/rhaensel/DynaTag/PDX_SCLC_DynaTag/Motifs_homer/denovo_results"
email="rhaensel@uni-koeln.de"
conda_env="/projects/ag-haensel/tools/.conda/envs/abc-model-env"
time_limit="24:00:00"
memory="32G"
cpus="8"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all normalized peak files
for normalized_bed in "${input_dir}"/*.bed; do
    # Extract base name for the job
    base_name=$(basename "${normalized_bed}" .bed)
    bg_file="${background_dir}/bg_kmer2_${base_name}.bg.positions.bed"
    motifs_output="${output_dir}/motifs_${base_name}"

    # Validate if background file exists before submitting the motif analysis
    if [[ ! -f "$bg_file" ]]; then
        echo "Error: Background file $bg_file not found for $base_name. Skipping."
        continue
    fi

    # Submit SLURM job for de novo motif discovery
    slurm_script="#!/bin/bash
#SBATCH --job-name=denovo_${base_name}
#SBATCH --output=${motifs_output}_slurm.out
#SBATCH --error=${motifs_output}_slurm.err
#SBATCH --time=${time_limit}
#SBATCH --mem=${memory}
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${email}

conda activate ${conda_env}

findMotifsGenome.pl \\
    ${normalized_bed} \\
    /projects/ag-haensel/Pascal/genome_files/hg38_bowtie/Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
    ${motifs_output} \\
    -size given \\
    -len 8,10,12 \\
    -S 25 \\
    -mis 2 \\
    -mset vertebrates \\
    -gc \\
    -bg ${bg_file} \\
    -redundant 2.0

conda deactivate
"
    slurm_script_file="submit_denovo_${base_name}.sbatch"
    echo "$slurm_script" > "$slurm_script_file"
    sbatch "$slurm_script_file"
done

echo "De novo motif discovery jobs submitted."
```
