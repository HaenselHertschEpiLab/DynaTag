# Response_16
## Response Fig. 2a Tornado plots of NANOG, SOX2 and OCT4 using matched ESC DynaTag, CUT&RUN and ChIP-seq data sets.
```bash
# Performed with Galaxy Europe server
# ChIP-seq NANOG
#!/bin/bash
# Galaxy computeMatrix script
# Ensure proper symbolic linking and computeMatrix execution

# Create symbolic links for input files
ln -f -s '/data/dnb10/galaxy_db/files/4/a/6/dataset_4a629cfc-d74f-4047-bb79-ea9392e73813.dat' \
  'SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean_cpm.bw_0.bw'

ln -f -s '/data/dnb10/galaxy_db/files/b/5/c/dataset_b5c1222c-e8c0-4bd1-b9c9-46144805f611.dat' \
  'SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw_1.bw'

ln -f -s '/data/dnb10/galaxy_db/files/4/9/2/dataset_492ac629-f4e0-405c-9471-0e6f8999ceec.dat' \
  'Nanog_mm10_target_genes_10kb_bed_0.bed'

# Execute computeMatrix
computeMatrix reference-point \
  --regionsFileName 'Nanog_mm10_target_genes_10kb_bed_0.bed' \
  --scoreFileName 'SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean_cpm.bw_0.bw' \
                 'SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw_1.bw' \
  --outFileName '/data/jwd02f/main/076/600/76600282/outputs/dataset_0f44ad38-30d1-475e-a123-dcf9733c9c4c.dat' \
  --samplesLabel SRR10992265_NANOG_ChIPseq SRR10992264_Input_ChIPseq \
  --numberOfProcessors "${GALAXY_SLOTS:-4}" \
  --referencePoint TSS \
  --beforeRegionStartLength 2500 \
  --afterRegionStartLength 2500 \
  --sortRegions 'keep' \
  --sortUsing 'mean' \
  --averageTypeBins 'mean' \
  --missingDataAsZero \
  --binSize 50

#!/bin/bash
# Galaxy plotHeatmap script
# Generate a heatmap using DeepTools plotHeatmap

plotHeatmap \
  --matrixFile '/data/dnb10/galaxy_db/files/0/f/4/dataset_0f44ad38-30d1-475e-a123-dcf9733c9c4c.dat' \
  --outFileName '/data/jwd05e/main/076/600/76600968/outputs/dataset_67a93446-fdd6-4f8b-8daf-9082122adf00.dat' \
  --plotFileFormat 'pdf' \
  --dpi '200' \
  --sortRegions 'no' \
  --sortUsing 'mean' \
  --averageTypeSummaryPlot 'mean' \
  --plotType 'lines' \
  --missingDataColor 'black' \
  --colorMap 'binary' \
  --alpha '1.0' \
  --zMin 0 \
  --xAxisLabel 'Distance' \
  --yAxisLabel 'Target Genes' \
  --heatmapWidth 7.5 \
  --heatmapHeight 7.5 \
  --whatToShow 'heatmap and colorbar' \
  --startLabel 'Upstream' \
  --endLabel 'Downstream' \
  --refPointLabel 'TSS' \
  --samplesLabel SRR10992265_NANOG_ChIPseq SRR10992264_Input_ChIPseq \
  --plotTitle 'NANOG ChIPSeq in NANOG ChIP-seq target genes plus.minus 2500bp' \
  --legendLocation 'best' \
  --labelRotation '0'

#!/bin/bash
# Galaxy computeMatrix and symbolic link script
# Create symbolic links and run computeMatrix

# Create symbolic links for input files
ln -f -s '/data/dnb10/galaxy_db/files/4/a/6/dataset_4a629cfc-d74f-4047-bb79-ea9392e73813.dat' \
  'SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean_cpm.bw_0.bw'

ln -f -s '/data/dnb10/galaxy_db/files/b/5/c/dataset_b5c1222c-e8c0-4bd1-b9c9-46144805f611.dat' \
  'SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw_1.bw'

ln -f -s '/data/dnb10/galaxy_db/files/4/9/2/dataset_492ac629-f4e0-405c-9471-0e6f8999ceec.dat' \
  'Nanog_mm10_target_genes_10kb_bed_0.bed'

# Run computeMatrix
computeMatrix reference-point \
  --regionsFileName 'Nanog_mm10_target_genes_10kb_bed_0.bed' \
  --scoreFileName 'SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean_cpm.bw_0.bw' \
                 'SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw_1.bw' \
  --outFileName '/data/jwd02f/main/076/600/76600974/outputs/dataset_90b84778-6a81-47f4-83a9-0644dfe08573.dat' \
  --samplesLabel SRR10992265_NANOG_ChIPseq SRR10992264_Input_ChIPseq \
  --numberOfProcessors "${GALAXY_SLOTS:-4}" \
  --referencePoint TSS \
  --beforeRegionStartLength 500 \
  --afterRegionStartLength 500 \
  --sortRegions 'keep' \
  --sortUsing 'mean' \
  --averageTypeBins 'mean' \
  --missingDataAsZero \
  --binSize 50

#!/bin/bash
# Galaxy plotHeatmap script
# Generate a heatmap using DeepTools plotHeatmap

plotHeatmap \
  --matrixFile '/data/dnb10/galaxy_db/files/9/0/b/dataset_90b84778-6a81-47f4-83a9-0644dfe08573.dat' \
  --outFileName '/data/jwd05e/main/076/601/76601004/outputs/dataset_acb8fb59-7e63-492e-ac3e-e998c64ab96c.dat' \
  --plotFileFormat 'pdf' \
  --dpi '200' \
  --sortRegions 'no' \
  --sortUsing 'mean' \
  --averageTypeSummaryPlot 'mean' \
  --plotType 'lines' \
  --missingDataColor 'black' \
  --colorMap 'binary' \
  --alpha '1.0' \
  --zMin 0 \
  --xAxisLabel 'Distance' \
  --yAxisLabel 'Target Genes' \
  --heatmapWidth 7.5 \
  --heatmapHeight 7.5 \
  --whatToShow 'heatmap and colorbar' \
  --startLabel 'Upstream' \
  --endLabel 'Downstream' \
  --refPointLabel 'TSS' \
  --samplesLabel SRR10992265_NANOG_ChIPseq SRR10992264_Input_ChIPseq \
  --plotTitle 'NANOG ChIPSeq, Input in NANOG ChIP-seq target genes plus.minus 500bp' \
  --legendLocation 'best' \
  --labelRotation '0'

# ChIP-seq SOX2 - same parameters as above, input files SRR10992267_SOX2_ChIPseq_30000000_mm10_norm_clean_cpm.bw,
# SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw, SOX2_mm10_target.genes_10kb.bed

# ChIP-seq OCT4 - same parameters as above, input files SRR10992266_OCT4_ChIPseq_30000000_mm10_norm_clean_cpm.bw,
# SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw, OCT4_mm10_target.genes_10kb.bed

# CUT&RUN NANOG - same parameters as above, input files SRR10992265_NANOG_ChIPseq_30000000_mm10_norm_clean_cpm.bw,
# SRR10992264_WT_input_ChIPseq_30000000_mm10_norm_clean_cpm.bw, Nanog_mm10_target.genes_10kb.bed

# CUT&RUN SOX2 - same parameters as above, input files GSM7019060_SOX2_CUTnRUN_2i.mm10_yeast.bw,
# GSM7019061_SOX2_CUTnRUN_SL.mm10_yeast.bw, SOX2_mm10_target.genes_10kb.bed

# CUT&RUN OCT4 - same parameters as above, input files GSM7019039_OCT4_CUTnRUN_2i.mm10_yeast.bw,
# GSM7019040_OCT4_CUTnRUN_SL.mm10_yeast.bw, OCT4_mm10_target.genes_10kb.bed

# DynaTag NANOG - same parameters as above, input files ESC-NANOG-G1.mm10.merged_cpm.bw,
# ESC-NANOG-G2.mm10.merged_cpm.bw, ESC-NANOG-S.mm10.merged_cpm.bw, Nanog_mm10_target.genes_10kb.bed

# DynaTag SOX2 - same parameters as above, input files ESC-SOX2-G1.mm10.merged_cpm.bw,
# ESC-SOX2-G2.mm10.merged_cpm.bw, ESC-SOX2-S.mm10.merged_cpm.bw SOX2_mm10_target.genes_10kb.bed

# DynaTag OCT4 - same parameters as above, input files ESC-OCT4-G1.mm10.merged_cpm.bw,
# ESC-OCT4-G2.mm10.merged_cpm.bw, ESC-OCT4-S.mm10.merged_cpm.bw, OCT4_mm10_target.genes_10kb.bed
```
## Response Fig. 2b HOMER2 known motif analyses NANOG, SOX2 and OCT4 using matched ESC DynaTag, CUT&RUN and ChIP-seq data sets.
```bash
cat submit_normalize_bg_motifs.sh
#!/bin/bash

# Parameters
genome_fa="/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa"
output_dir="motifs_results"
background_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/bg_results"
email="rhaensel@uni-koeln.de"
conda_env="/projects/ag-haensel/tools/.conda/envs/abc-model-env"
kmer_value=2  # Set kmer to 2 for all jobs

# Create output directories
mkdir -p ${output_dir}
mkdir -p ${background_dir}

# Process all .bed files in the current directory
for input_bed in *.bed; do
    echo "Processing file: ${input_bed}"

    # Calculate average peak size
    average_peak_size=$(awk '{sum += $3 - $2} END {if (NR > 0) print int(sum / NR); else print 200}' ${input_bed})
    echo "Calculated average peak size for ${input_bed}: ${average_peak_size} bp"

    # Use the calculated average peak size as the fixed_size
    fixed_size=${average_peak_size}

    # Normalize BED file locally
    normalized_bed="normalized_${input_bed}"
    awk -v size=${fixed_size} 'BEGIN {OFS="\t"} {
        mid = int(($2 + $3) / 2);
        start = mid - int(size / 2);
        end = mid + int(size / 2);
        if (start < 0) start = 0;
        print $1, start, end, ".", "1000";
    }' ${input_bed} > ${normalized_bed}

    echo "Step 1: Normalization completed for ${input_bed}."

    # Step 2: Generate Background Files
    bg_output_prefix="${background_dir}/bg_kmer${kmer_value}_${input_bed%.bed}"

    slurm_script="#!/bin/bash
#SBATCH --job-name=bg_kmer${kmer_value}_${input_bed%.bed}
#SBATCH --output=${bg_output_prefix}_slurm.out
#SBATCH --error=${bg_output_prefix}_slurm.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${email}

conda activate ${conda_env}

homer2 background \\
    -p ${normalized_bed} \\
    -g ${genome_fa} \\
    -pkmer ${kmer_value} \\
    -N 100000 \\
    -NN 100000000 \\
    -allowTargetOverlaps \\
    -allowBgOverlaps \\
    -o ${bg_output_prefix}

conda deactivate
"
    slurm_script_file="submit_bg_kmer${kmer_value}_${input_bed%.bed}.sbatch"
    echo "$slurm_script" > "$slurm_script_file"
    bg_job_id=$(sbatch --parsable "$slurm_script_file")
    echo "Step 2: Background file generation submitted for ${input_bed} with job ID: $bg_job_id"

    # Step 3: Submit Motif Analysis Job
    bg_file="${background_dir}/bg_kmer${kmer_value}_${input_bed%.bed}.bg.positions.bed"
    motifs_output="${output_dir}/motifs_kmer${kmer_value}_${input_bed%.bed}"

    # Validate if background file exists before submitting the motif analysis
    if [[ ! -f "$bg_file" ]]; then
        echo "Error: Background file $bg_file not found for kmer ${kmer_value}. Exiting."
        exit 1
    fi

    slurm_script="#!/bin/bash
#SBATCH --job-name=motifs_kmer${kmer_value}_${input_bed%.bed}
#SBATCH --output=${motifs_output}_slurm.out
#SBATCH --error=${motifs_output}_slurm.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${email}
#SBATCH --dependency=afterok:${bg_job_id}

conda activate ${conda_env}

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
    -redundant 2.0 \\
    -maxN 0.7 \\
    -nlen 3 \\
    -nmax 160 \\
    -e 0.0 \\
    -minlp -10.0

conda deactivate
"
    slurm_script_file="submit_motifs_kmer${kmer_value}_${input_bed%.bed}.sbatch"
    echo "$slurm_script" > "$slurm_script_file"
    sbatch "$slurm_script_file"

    echo "Step 3: Motif analysis job submitted for ${input_bed}."
done

cat normalize_bed_files.sh
#!/bin/bash

# Parameters
input_files="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/ESC-ATAC_master_peaks.mm10.bed /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/*_master_peaks.bed /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/GSM4291125_mESC_ChIPseq_peaks/*_peaks.bed /scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/peaks/peaks_GSE224292_mESC_CUTnRUN/*_peaks.bed"  # Process all BED files in the directory
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

cat generate_background.sh
#!/bin/bash

# Parameters
normalized_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/bg_results"
genome_fa="/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa"
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

cat run_motif_analysis.sh
#!/bin/bash

# Parameters
normalized_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/bg_results"
output_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/motifs_results"
genome_fa="/projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa"
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

cat run_denovo_motif_analysis.sh
#!/bin/bash

# Parameters
input_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/normalized_peaks"
background_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/bg_results"
output_dir="/scratch/rhaensel/DynaTag/ESC_EpiLC_DynaTag/Motifs_homer/denovo_results"
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
    /projects/ag-haensel/Pascal/genome_files/mm10_bowtie/mm10.fa \\
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


