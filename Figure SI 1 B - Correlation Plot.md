# Preprocessing in DynaTag_bulk_ESC_EpiLC.md
# plotCorrelation
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
