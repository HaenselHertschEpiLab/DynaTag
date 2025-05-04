# DiffBind Rscript
```R
SessionInfo

R version 3.5.1 (2018-07-02)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Rocky Linux 9.3 (Blue Onyx)

Matrix products: default
BLAS/LAPACK: /usr/local/tools/_conda/envs/mulled-v1-655069e4a80d4f19c337748eeb95d6ab9de06ebdab00b06603f42c6f080aae0b/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] rjson_0.2.20                DiffBind_2.10.0
 [3] SummarizedExperiment_1.12.0 DelayedArray_0.8.0
 [5] BiocParallel_1.16.6         matrixStats_0.54.0
 [7] Biobase_2.42.0              GenomicRanges_1.34.0
 [9] GenomeInfoDb_1.18.1         IRanges_2.16.0
[11] S4Vectors_0.20.1            BiocGenerics_0.28.0
[13] getopt_1.20.3

suppressPackageStartupMessages({
    library('getopt')
    library('DiffBind')
    library('rjson')
})

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'infile' , 'i', 1, "character",
    'outfile' , 'o', 1, "character",
    'scorecol', 'n', 1, "integer",
    'lowerbetter', 'l', 1, "logical",
    'summits', 's', 1, "integer",
    'th', 't', 1, "double",
    'format', 'f', 1, "character",
    'plots' , 'p', 2, "character",
    'bmatrix', 'b', 0, "logical",
    "rdaOpt", "r", 0, "logical",
    'infoOpt' , 'a', 0, "logical",
    'verbose', 'v', 2, "integer",
    'help' , 'h', 0, "logical"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

parser <- newJSONParser()
parser$addData(opt$infile)
factorList <- parser$getObject()
filenamesIn <- unname(unlist(factorList[[1]][[2]]))
peaks <- filenamesIn[grepl("peaks.bed", filenamesIn)]
bams <- filenamesIn[grepl("bamreads.bam", filenamesIn)]
ctrls <- filenamesIn[grepl("bamcontrol.bam", filenamesIn)]

# get the group and sample id from the peaks filenames
groups <- sapply(strsplit(peaks,"-"), `[`, 1)
samples <- sapply(strsplit(peaks,"-"), `[`, 2)

if ( length(ctrls) != 0 ) {
    sampleTable <- data.frame(SampleID=samples,
                        Condition=groups,
                        bamReads=bams,
                        bamControl=ctrls,
                        Peaks=peaks,
                        Tissue=samples) # using "Tissue" column to display ids as labels in PCA plot
} else {

    sampleTable <- data.frame(SampleID=samples,
                        Replicate=samples,
                        Condition=groups,
                        bamReads=bams,
                        Peaks=peaks,
                        Tissue=samples)
}

sample = dba(sampleSheet=sampleTable, peakFormat='bed', scoreCol=opt$scorecol, bLowerScoreBetter=opt$lowerbetter)

if ( !is.null(opt$summits) ) {
    sample_count = dba.count(sample, summits=opt$summits)
} else {
    sample_count = dba.count(sample)
}

sample_contrast = dba.contrast(sample_count, categories=DBA_CONDITION, minMembers=2)
sample_analyze = dba.analyze(sample_contrast)
diff_bind = dba.report(sample_analyze, th=opt$th)

# Generate plots
if ( !is.null(opt$plots) ) {
    pdf(opt$plots)
    orvals = dba.plotHeatmap(sample_analyze, contrast=1, correlations=FALSE, cexCol=0.8, th=opt$th)
    dba.plotPCA(sample_analyze, contrast=1, th=opt$th, label=DBA_TISSUE, labelSize=0.3)
    dba.plotMA(sample_analyze, th=opt$th)
    dba.plotVolcano(sample_analyze, th=opt$th)
    dba.plotBox(sample_analyze, th=opt$th)
    dev.off()
}

# Output differential binding sites
resSorted <- diff_bind[order(diff_bind$FDR),]
# Convert from GRanges (1-based) to 0-based format (adapted from https://www.biostars.org/p/89341/)
if (opt$format == "bed") {
    resSorted  <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)))
} else if (opt$format == "interval") {
     # Output as interval
    df <- as.data.frame(resSorted)
    extrainfo <- NULL
    for (i in 1:nrow(df)) {
        extrainfo[i] <- paste0(c(df$width[i], df[i, 6:ncol(df)]), collapse="|")
    }
    resSorted  <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)),
        Comment=extrainfo)
} else {
    # Output as 0-based tabular
    resSorted <- data.frame(Chrom=seqnames(resSorted),
        Start=start(resSorted) - 1,
        End=end(resSorted),
        Name=rep("DiffBind", length(resSorted)),
        Score=rep("0", length(resSorted)),
        Strand=gsub("\\*", ".", strand(resSorted)),
        mcols(resSorted))
}
write.table(resSorted, file = opt$outfile, sep="\t", quote = FALSE, row.names = FALSE)

# Output binding affinity scores
if (!is.null(opt$bmatrix)) {
    bmat <- dba.peakset(sample_count, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
    # Output as 0-based tabular
    bmat <- data.frame(Chrom=bmat[, 1],
        Start=bmat[, 2] - 1,
        End=bmat[, 3],
        bmat[, 4:ncol(bmat)])
    write.table(bmat, file="bmatrix.tab", sep="\t", quote=FALSE, row.names=FALSE)
}

# Output RData file
if (!is.null(opt$rdaOpt)) {
    save.image(file = "DiffBind_analysis.RData")
}

# Output analysis info
if (!is.null(opt$infoOpt)) {
    info <- "DiffBind_analysis_info.txt"
    cat("dba.count Info\n\n", file=info, append = TRUE)
    capture.output(sample, file=info, append=TRUE)
    cat("\ndba.analyze Info\n\n", file=info, append = TRUE)
    capture.output(sample_analyze, file=info, append=TRUE)
    cat("\nSessionInfo\n\n", file=info, append = TRUE)
    capture.output(sessionInfo(), file=info, append=TRUE)
}

#DiffBind_session_info
#FOXA1
dba.count Info

4 Samples, 59355 sites in matrix (318214 total):
                  ID             Tissue          Condition Replicate Caller
1 FOXA1_CHEM_2_peaks FOXA1_CHEM_2_peaks FOXA1_Chemotherapy         2    raw
2 FOXA1_CHEM_1_peaks FOXA1_CHEM_1_peaks FOXA1_Chemotherapy         1    raw
3 FOXA1_CTRL_2_peaks FOXA1_CTRL_2_peaks      FOXA1_Control         4    raw
4 FOXA1_CTRL_1_peaks FOXA1_CTRL_1_peaks      FOXA1_Control         3    raw
  Intervals
1    100656
2     98775
3    116643
4    117450

dba.analyze Info

4 Samples, 59355 sites in matrix:
                  ID             Tissue          Condition Replicate Caller
1 FOXA1_CHEM_2_peaks FOXA1_CHEM_2_peaks FOXA1_Chemotherapy         2 counts
2 FOXA1_CHEM_1_peaks FOXA1_CHEM_1_peaks FOXA1_Chemotherapy         1 counts
3 FOXA1_CTRL_2_peaks FOXA1_CTRL_2_peaks      FOXA1_Control         4 counts
4 FOXA1_CTRL_1_peaks FOXA1_CTRL_1_peaks      FOXA1_Control         3 counts
  Intervals FRiP
1     59355 0.31
2     59355 0.32
3     59355 0.37
4     59355 0.36

1 Contrast:
              Group1 Members1        Group2 Members2 DB.DESeq2
1 FOXA1_Chemotherapy        2 FOXA1_Control        2      1874



#ASCL1
dba.count Info

4 Samples, 40536 sites in matrix (110950 total):
                  ID             Tissue          Condition Replicate Caller
1 ASCL1_CHEM_2_peaks ASCL1_CHEM_2_peaks ASCL1_Chemotherapy         2    raw
2 ASCL1_CHEM_1_peaks ASCL1_CHEM_1_peaks ASCL1_Chemotherapy         1    raw
3 ASCL1_CTRL_2_peaks ASCL1_CTRL_2_peaks      ASCL1_Control         4    raw
4 ASCL1_CTRL_1_peaks ASCL1_CTRL_1_peaks      ASCL1_Control         3    raw
  Intervals
1     51588
2     51021
3     50358
4     48906

dba.analyze Info

4 Samples, 40536 sites in matrix:
                  ID             Tissue          Condition Replicate Caller
1 ASCL1_CHEM_2_peaks ASCL1_CHEM_2_peaks ASCL1_Chemotherapy         2 counts
2 ASCL1_CHEM_1_peaks ASCL1_CHEM_1_peaks ASCL1_Chemotherapy         1 counts
3 ASCL1_CTRL_2_peaks ASCL1_CTRL_2_peaks      ASCL1_Control         4 counts
4 ASCL1_CTRL_1_peaks ASCL1_CTRL_1_peaks      ASCL1_Control         3 counts
  Intervals FRiP
1     40536 0.29
2     40536 0.29
3     40536 0.42
4     40536 0.42

1 Contrast:
              Group1 Members1        Group2 Members2 DB.DESeq2
1 ASCL1_Chemotherapy        2 ASCL1_Control        2      7301



MYC
dba.count Info

4 Samples, 56930 sites in matrix (277877 total):
                ID           Tissue        Condition Replicate Caller Intervals
1 MYC_CHEM_2_peaks MYC_CHEM_2_peaks MYC_Chemotherapy         2    raw    115321
2 MYC_CHEM_1_peaks MYC_CHEM_1_peaks MYC_Chemotherapy         1    raw    111100
3 MYC_CTRL_2_peaks MYC_CTRL_2_peaks      MYC_Control         4    raw     88603
4 MYC_CTRL_1_peaks MYC_CTRL_1_peaks      MYC_Control         3    raw     79305

dba.analyze Info

4 Samples, 56930 sites in matrix:
                ID           Tissue        Condition Replicate Caller Intervals
1 MYC_CHEM_2_peaks MYC_CHEM_2_peaks MYC_Chemotherapy         2 counts     56930
2 MYC_CHEM_1_peaks MYC_CHEM_1_peaks MYC_Chemotherapy         1 counts     56930
3 MYC_CTRL_2_peaks MYC_CTRL_2_peaks      MYC_Control         4 counts     56930
4 MYC_CTRL_1_peaks MYC_CTRL_1_peaks      MYC_Control         3 counts     56930
  FRiP
1 0.34
2 0.34
3 0.34
4 0.36

1 Contrast:
            Group1 Members1      Group2 Members2 DB.DESeq2
1 MYC_Chemotherapy        2 MYC_Control        2      1281

#NEUROD1
dba.count Info

4 Samples, 57114 sites in matrix (451612 total):
                    ID               Tissue            Condition Replicate
1 NEUROD1_CHEM_2_peaks NEUROD1_CHEM_2_peaks NEUROD1_Chemotherapy         2
2 NEUROD1_CHEM_1_peaks NEUROD1_CHEM_1_peaks NEUROD1_Chemotherapy         1
3 NEUROD1_CTRL_2_peaks NEUROD1_CTRL_2_peaks      NEUROD1_Control         4
4 NEUROD1_CTRL_1_peaks NEUROD1_CTRL_1_peaks      NEUROD1_Control         3
  Caller Intervals
1    raw    147759
2    raw    147024
3    raw    125704
4    raw    124406

dba.analyze Info

4 Samples, 57114 sites in matrix:
                    ID               Tissue            Condition Replicate
1 NEUROD1_CHEM_2_peaks NEUROD1_CHEM_2_peaks NEUROD1_Chemotherapy         2
2 NEUROD1_CHEM_1_peaks NEUROD1_CHEM_1_peaks NEUROD1_Chemotherapy         1
3 NEUROD1_CTRL_2_peaks NEUROD1_CTRL_2_peaks      NEUROD1_Control         4
4 NEUROD1_CTRL_1_peaks NEUROD1_CTRL_1_peaks      NEUROD1_Control         3
  Caller Intervals FRiP
1 counts     57114 0.33
2 counts     57114 0.34
3 counts     57114 0.41
4 counts     57114 0.41

1 Contrast:
                Group1 Members1          Group2 Members2 DB.DESeq2
1 NEUROD1_Chemotherapy        2 NEUROD1_Control        2      1117

dba.count Info

4 Samples, 52959 sites in matrix (210483 total):
                 ID            Tissue         Condition Replicate Caller
1 NRF1_CHEM_2_peaks NRF1_CHEM_2_peaks NRF1_Chemotherapy         2    raw
2 NRF1_CHEM_1_peaks NRF1_CHEM_1_peaks NRF1_Chemotherapy         1    raw
3 NRF1_CTRL_2_peaks NRF1_CTRL_2_peaks      NRF1_Control         4    raw
4 NRF1_CTRL_1_peaks NRF1_CTRL_1_peaks      NRF1_Control         3    raw
  Intervals
1     84184
2     81623
3     79976
4     80673

#NRF1
dba.analyze Info

4 Samples, 52959 sites in matrix:
                 ID            Tissue         Condition Replicate Caller
1 NRF1_CHEM_2_peaks NRF1_CHEM_2_peaks NRF1_Chemotherapy         2 counts
2 NRF1_CHEM_1_peaks NRF1_CHEM_1_peaks NRF1_Chemotherapy         1 counts
3 NRF1_CTRL_2_peaks NRF1_CTRL_2_peaks      NRF1_Control         4 counts
4 NRF1_CTRL_1_peaks NRF1_CTRL_1_peaks      NRF1_Control         3 counts
  Intervals FRiP
1     52959 0.37
2     52959 0.38
3     52959 0.47
4     52959 0.46

1 Contrast:
             Group1 Members1       Group2 Members2 DB.DESeq2
1 NRF1_Chemotherapy        2 NRF1_Control        2      1733


#POU2F3
dba.count Info

4 Samples, 9891 sites in matrix (135086 total):
                   ID              Tissue           Condition Replicate Caller
1 POU2F3_CHEM_2_peaks POU2F3_CHEM_2_peaks POU2F3_Chemotherapy         2    raw
2 POU2F3_CHEM_1_peaks POU2F3_CHEM_1_peaks POU2F3_Chemotherapy         1    raw
3 POU2F3_CTRL_2_peaks POU2F3_CTRL_2_peaks      POU2F3_Control         4    raw
4 POU2F3_CTRL_1_peaks POU2F3_CTRL_1_peaks      POU2F3_Control         3    raw
  Intervals
1      3127
2      3151
3    131526
4     11032

dba.analyze Info

4 Samples, 9891 sites in matrix:
                   ID              Tissue           Condition Replicate Caller
1 POU2F3_CHEM_2_peaks POU2F3_CHEM_2_peaks POU2F3_Chemotherapy         2 counts
2 POU2F3_CHEM_1_peaks POU2F3_CHEM_1_peaks POU2F3_Chemotherapy         1 counts
3 POU2F3_CTRL_2_peaks POU2F3_CTRL_2_peaks      POU2F3_Control         4 counts
4 POU2F3_CTRL_1_peaks POU2F3_CTRL_1_peaks      POU2F3_Control         3 counts
  Intervals FRiP
1      9891 0.13
2      9891 0.13
3      9891 0.24
4      9891 0.24

1 Contrast:
               Group1 Members1         Group2 Members2 DB.DESeq2
1 POU2F3_Chemotherapy        2 POU2F3_Control        2      3999

#NFIB
dba.count Info

4 Samples, 46335 sites in matrix (100318 total):
                 ID            Tissue         Condition Replicate Caller
1 NFIB_CHEM_2_peaks NFIB_CHEM_2_peaks NFIB_Chemotherapy         2    raw
2 NFIB_CHEM_1_peaks NFIB_CHEM_1_peaks NFIB_Chemotherapy         1    raw
3 NFIB_CTRL_2_peaks NFIB_CTRL_2_peaks      NFIB_Control         4    raw
4 NFIB_CTRL_1_peaks NFIB_CTRL_1_peaks      NFIB_Control         3    raw
  Intervals
1     51918
2     52134
3     52982
4     52593

dba.analyze Info

4 Samples, 46335 sites in matrix:
                 ID            Tissue         Condition Replicate Caller
1 NFIB_CHEM_2_peaks NFIB_CHEM_2_peaks NFIB_Chemotherapy         2 counts
2 NFIB_CHEM_1_peaks NFIB_CHEM_1_peaks NFIB_Chemotherapy         1 counts
3 NFIB_CTRL_2_peaks NFIB_CTRL_2_peaks      NFIB_Control         4 counts
4 NFIB_CTRL_1_peaks NFIB_CTRL_1_peaks      NFIB_Control         3 counts
  Intervals FRiP
1     46335 0.37
2     46335 0.37
3     46335 0.44
4     46335 0.44

1 Contrast:
             Group1 Members1       Group2 Members2 DB.DESeq2
1 NFIB_Chemotherapy        2 NFIB_Control        2      2081

#p53
dba.count Info

4 Samples, 47823 sites in matrix (389275 total):
                ID           Tissue        Condition Replicate Caller Intervals
1 p53_CHEM_2_peaks p53_CHEM_2_peaks p53_Chemotherapy         2    raw    125089
2 p53_CHEM_1_peaks p53_CHEM_1_peaks p53_Chemotherapy         1    raw    120002
3 p53_CTRL_2_peaks p53_CTRL_2_peaks      p53_Control         4    raw    114455
4 p53_CTRL_1_peaks p53_CTRL_1_peaks      p53_Control         3    raw    116373

dba.analyze Info

4 Samples, 47823 sites in matrix:
                ID           Tissue        Condition Replicate Caller Intervals
1 p53_CHEM_2_peaks p53_CHEM_2_peaks p53_Chemotherapy         2 counts     47823
2 p53_CHEM_1_peaks p53_CHEM_1_peaks p53_Chemotherapy         1 counts     47823
3 p53_CTRL_2_peaks p53_CTRL_2_peaks      p53_Control         4 counts     47823
4 p53_CTRL_1_peaks p53_CTRL_1_peaks      p53_Control         3 counts     47823
  FRiP
1 0.25
2 0.25
3 0.32
4 0.32

1 Contrast:
            Group1 Members1      Group2 Members2 DB.DESeq2
1 p53_Chemotherapy        2 p53_Control        2       575

#YAP1
dba.count Info

4 Samples, 46861 sites in matrix (475618 total):
                 ID            Tissue         Condition Replicate Caller
1 YAP1_CHEM_2_peaks YAP1_CHEM_2_peaks YAP1_Chemotherapy         2    raw
2 YAP1_CHEM_1_peaks YAP1_CHEM_1_peaks YAP1_Chemotherapy         1    raw
3 YAP1_CTRL_2_peaks YAP1_CTRL_2_peaks      YAP1_Control         4    raw
4 YAP1_CTRL_1_peaks YAP1_CTRL_1_peaks      YAP1_Control         3    raw
  Intervals
1    156579
2    156471
3    114073
4    115126

dba.analyze Info

4 Samples, 46861 sites in matrix:
                 ID            Tissue         Condition Replicate Caller
1 YAP1_CHEM_2_peaks YAP1_CHEM_2_peaks YAP1_Chemotherapy         2 counts
2 YAP1_CHEM_1_peaks YAP1_CHEM_1_peaks YAP1_Chemotherapy         1 counts
3 YAP1_CTRL_2_peaks YAP1_CTRL_2_peaks      YAP1_Control         4 counts
4 YAP1_CTRL_1_peaks YAP1_CTRL_1_peaks      YAP1_Control         3 counts
  Intervals FRiP
1     46861 0.28
2     46861 0.28
3     46861 0.36
4     46861 0.36

1 Contrast:
             Group1 Members1       Group2 Members2 DB.DESeq2
1 YAP1_Chemotherapy        2 YAP1_Control        2       638


#p73
dba.count Info

4 Samples, 29430 sites in matrix (130498 total):
                ID           Tissue        Condition Replicate Caller Intervals
1 p73_CHEM_2_peaks p73_CHEM_2_peaks p73_Chemotherapy         2    raw     54857
2 p73_CHEM_1_peaks p73_CHEM_1_peaks p73_Chemotherapy         1    raw     54076
3 p73_CTRL_2_peaks p73_CTRL_2_peaks      p73_Control         4    raw     40854
4 p73_CTRL_1_peaks p73_CTRL_1_peaks      p73_Control         3    raw     40607

dba.analyze Info

4 Samples, 29430 sites in matrix:
                ID           Tissue        Condition Replicate Caller Intervals
1 p73_CHEM_2_peaks p73_CHEM_2_peaks p73_Chemotherapy         2 counts     29430
2 p73_CHEM_1_peaks p73_CHEM_1_peaks p73_Chemotherapy         1 counts     29430
3 p73_CTRL_2_peaks p73_CTRL_2_peaks      p73_Control         4 counts     29430
4 p73_CTRL_1_peaks p73_CTRL_1_peaks      p73_Control         3 counts     29430
  FRiP
1 0.33
2 0.33
3 0.49
4 0.49

1 Contrast:
            Group1 Members1      Group2 Members2 DB.DESeq2
1 p73_Chemotherapy        2 p73_Control        2      1435


#SP2
dba.count Info

4 Samples, 36314 sites in matrix (417668 total):
                ID           Tissue        Condition Replicate Caller Intervals
1 SP2_CHEM_2_peaks SP2_CHEM_2_peaks SP2_Chemotherapy         2    raw    158850
2 SP2_CHEM_1_peaks SP2_CHEM_1_peaks SP2_Chemotherapy         1    raw    157277
3 SP2_CTRL_2_peaks SP2_CTRL_2_peaks      SP2_Control         4    raw     77469
4 SP2_CTRL_1_peaks SP2_CTRL_1_peaks      SP2_Control         3    raw     78497

dba.analyze Info

4 Samples, 36314 sites in matrix:
                ID           Tissue        Condition Replicate Caller Intervals
1 SP2_CHEM_2_peaks SP2_CHEM_2_peaks SP2_Chemotherapy         2 counts     36314
2 SP2_CHEM_1_peaks SP2_CHEM_1_peaks SP2_Chemotherapy         1 counts     36314
3 SP2_CTRL_2_peaks SP2_CTRL_2_peaks      SP2_Control         4 counts     36314
4 SP2_CTRL_1_peaks SP2_CTRL_1_peaks      SP2_Control         3 counts     36314
  FRiP
1 0.28
2 0.28
3 0.43
4 0.43

1 Contrast:
            Group1 Members1      Group2 Members2 DB.DESeq2
1 SP2_Chemotherapy        2 SP2_Control        2       543
```

# Filter DBA tables with awk to retrieve increased and decreased occupancy
```bash
for file in *.txt; do   echo "Processing $file";   awk 'NR>1 {
    if ($1 != "MT") {
      $1 = "chr" $1;
      $12 = -log($12) / log(10);
      if ($12 > 1.3 && $10 < -0.5) {
        print
      }
    }
  }' "$file" | awk 'BEGIN {OFS="\t"} {$1=$1; print}' | sort -k1,1 -k2,2n > "${file%_FDR0.1.txt}_FDR.below.0.05_FC.below.0.5_Chemo.lost_sort.bed"; done

  for file in *.txt; do   echo "Processing $file";   awk 'NR>1 {
    if ($1 != "MT") {
      $1 = "chr" $1;
      $12 = -log($12) / log(10);
      if ($12 > 1.3 && $10 > 0.5) {
        print
      }
    }
  }' "$file" | awk 'BEGIN {OFS="\t"} {$1=$1; print}' | sort -k1,1 -k2,2n > "${file%_FDR0.1.txt}_FDR.below.0.05_FC.above.0.5_Chemo.gained_sort.bed"; done
```
