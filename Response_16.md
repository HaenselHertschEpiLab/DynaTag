# Response_16
## Response Fig. 2a Tornado plots of NANOG, SOX2 and OCT4 uisng matched ESC DynaTag, CUT&RUN and ChIP-seq dat sets.
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







