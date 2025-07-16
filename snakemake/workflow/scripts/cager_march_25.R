# Load libraries
suppressPackageStartupMessages({
  library(CAGEfightR)
  library(CAGEr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(ggplot2)
  library(rtracklayer)
  library(BSgenome)
  library(BSgenome.Taestivum.ChineseSpring)
  library(glue)
  library(dplyr)
})

setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data")
output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps"
input_files <- c("SRO1.CS.se.star.bam", "SRO2.CS.se.star.bam", "SRO3.CS.se.star.bam",
                "SLE1.CS.se.star.bam", "SLE2.CS.se.star.bam", "SLE3.CS.se.star.bam",
                "SSP1.CS.se.star.bam", "SSP2.CS.se.star.bam", "SSP3.CS.se.star.bam",
                "SIS1.CS.se.star.bam", "SIS2.CS.se.star.bam", "SIS3.CS.se.star.bam",
                "CR1.CS.se.star.bam", "CR2.CS.se.star.bam", "CR3.CS.se.star.bam", "CR4.CS.se.star.bam", "CR5.CS.se.star.bam",
                "CL1.CS.se.star.bam", "CL2.CS.se.star.bam", "CL3.CS.se.star.bam", "CL4.CS.se.star.bam", "CL5.CS.se.star.bam")
sampleLabels <- c("SRO1", "SRO2", "SRO3", "SLE1", "SLE2", "SLE3", "SSP1", "SSP2", "SSP3", "SIS1", "SIS2", "SIS3", "CR1", "CR2", "CR3", "CR4", "CR5", "CL1", "CL2", "CL3", "CL4", "CL5")
genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                     inputFiles = input_files, 
                     inputFilesType = "bam",
                     sampleLabels = sampleLabels)
print(myCAGEset)

#myCAGEset <- getCTSS(myCAGEset, correctSystematicG = TRUE)
# Save myCAGEset object as it takes a long time to generate
#saveRDS(myCAGEset, "myCAGEset_correctsystematicG.rds")
# Load myCAGEset object
myCAGEset <- readRDS("myCAGEset_correctsystematicG.rds")
print(myCAGEset)
# Normalize to TPM
myCAGEset <- normalizeTagCount(myCAGEset, method =  "simpleTpm")
# filter - optional step
#       #myCAGEset <- filterLowExpCTSS(myCAGEset, thresholdIsTpm = TRUE, threshold = 0.5) # Flag CTSSs with signal < 0.5 TPM - nrPassThreshold in future
#       #filtered <- filteredCTSSidx(myCAGEset)
#       #sum(filtered)  # Count of CTSSs retained - [1] 898367
#       #sum(!filtered)  # Count of CTSSs removed - [1] 12103278
#       #filtered_idx <- which(filtered)  # Convert to numeric index
#       #myCAGEset <- myCAGEset[filtered_idx, ]
##############

# Cluster CTSS into transcriptional clusters (TCs)
# distclu
distclu_myCAGEset <- distclu(myCAGEset, maxDist = 50, keepSingletonsAbove = 5)
saveRDS(distclu_myCAGEset, "no_merge_reps/no_merge_distclu_md20_single_5_myCAGEset.rds")
# Read in distclu object

#distclu_myCAGEset <- readRDS("no_merge_reps/distclu_myCAGEset.rds")

#convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
#  # This function is required for the annotateCTSS function to work!
#  print('Generating own granges object from GFF file .... ')
#  data_frame <- read.delim(gff_file, header=F, skip = 1)
#  data_frame <- data_frame %>%
#    mutate(extracted_id = ifelse(
#    grepl('=([^\\.]+)\\.', V9),
#    sub('.*=([^\\.]+)\\..*', '\\1', V9),
#    NA)
#  )
#  # Replace "mRNA" with "transcript" - Thiis is important for the plot annot to work!
#  data_frame$V3[data_frame$V3 == "mRNA"] <- "transcript"
#  seqnames <- Rle(data_frame$V1) # e.g. Chr1A 
#  targetRanges <- IRanges(data_frame$V4, data_frame$V5) 
#  strands <- Rle(data_frame$V7)
#  gene_name <- Rle(data_frame$extracted_id)
#  type <- Rle(data_frame$V3)
#  
#  transcript_type <- Rle(rep(transcript_type, length(seqnames)))
#  gr <- GRanges(seqnames = seqnames, ranges = targetRanges, strand = strands, type = type, transcript_type = transcript_type, gene_name = gene_name)
#  return(gr)
#}


#convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
#  print('Generating own GRanges object from GFF file .... ')
#  data_frame <- read.delim(gff_file, header = FALSE, skip = 1, stringsAsFactors = FALSE, colClasses = c(rep("character", 8), "character"))
#  colnames(data_frame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
#  data_frame$feature[data_frame$feature == "mRNA"] <- "transcript"
#  if (!is.character(data_frame$attributes)) {
#    data_frame$attributes <- as.character(data_frame$attributes)
#  }
#  data_frame <- data_frame %>%
#    mutate(
#      gene_name = case_when(
#        grepl("ID=([^;]+)", attributes) ~ sub(".*ID=([^;]+).*", "\\1", attributes),
#        grepl("gene_id=([^;]+)", attributes) ~ sub(".*gene_id=([^;]+).*", "\\1", attributes),
#        grepl("Parent=([^;]+)", attributes) ~ sub(".*Parent=([^;]+).*", "\\1", attributes),
#        TRUE ~ NA_character_  # Assign NA if no match is found
#      )
#    )
#  gr <- GRanges(
#    seqnames = Rle(data_frame$seqname),
#    ranges = IRanges(start = as.numeric(data_frame$start), end = as.numeric(data_frame$end)),
#    strand = Rle(data_frame$strand),
#    type = Rle(data_frame$feature),
#    transcript_type = Rle(rep(transcript_type, nrow(data_frame))),
#    gene_name = Rle(data_frame$gene_name)
#  )
#}
#granges <- convert_gff_to_granges("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")
## Check that granages contain some key genes in gene_name column 
## TraesCS7D03G1294200, TraesCS5B03G0371600.3, TraesCS5B03G0193200.1
#print('Debugging the granges object by checkiing it has these three genes :TraesCS7D03G1294200 TraesCS5B03G0371600.3 TraesCS5B03G0193200.1 TraesCS7D03G0742800.1 ')
#print(granges[granges$gene_name %in% c("TraesCS7D03G1294200", "TraesCS5B03G0371600.3", "TraesCS5B03G0193200.1", "TraesCS7D03G0742800.1", "TraesCS1A03G0000600", "TraesCS1A03G0000600.1")])
#distclu_myCAGEset <- annotateCTSS(distclu_myCAGEset, granges, upstream=500, downstream=500)
## Check if there are any occurances of "" as the gene name in the granges object
#print('Checking if there are any occurances of "" as the gene name in the granges object')
#print(granges[granges$gene_name == ""])
#print(print(table(granges$type)))
#print('Trying, remove genes so that just transcripts are used in annotation')
#granges <- granges[granges$type != "gene"]
#
#
##print(glue("Annotated CTSS data frame:"))
##print(colData(distclu_myCAGEset)[,c("librarySizes", "promoter", "exon", "intron", "unknown")])
##plot_file_path <- file.path('no_merge_reps', "CTSS_genomic_locations_500_500.pdf")
##plot <- plotAnnot(
##    distclu_myCAGEset,
##    scope = "counts",
##    title = "CTSS Genomic Location")
##ggsave(plot_file_path, plot, width = 8, height = 6, dpi = 900)
##message(glue("Saved the annotation plot to {plot_file_path}")) 
#
## Save all Annotated TagClusters to a csv
#annotate <- annotateTagClusters(distclu_myCAGEset, granges, upstream = 500, downstream = 500)
#for (sample in sampleLabels) {
#    message(glue("{Sys.Date()} - "))
#    sample_tagclusters <- tagClustersGR(annotate, sample = sample)
#    missing_gene_promoters <- sample_tagclusters[sample_tagclusters$annotation == "promoter" & sample_tagclusters$genes == ""]
#    # Print summary
#    print(glue("Total promoters without transcript ID: {length(missing_gene_promoters)}"))
#    # Debug statements 
#    tc_region <- GRanges(seqnames = "Chr1A", ranges = IRanges(start = 95585, end = 95686), strand = "+")
#    tc_annotation <- ranges2annot(tc_region, granges, upstream = 500, downstream = 500)
#    print(glue("Annotation assigned: {tc_annotation}"))
#    gene_name <- ranges2genes_extended(tc_region, granges, upstream = 500, downstream = 500)
#    print(glue("Gene name assigned: {gene_name}"))
#    # Debug finished 
#    exit()
#    if (length(overlapping_genes) > 0) {
#      print(granges[subjectHits(overlapping_genes)])
#    } else {
#        print("No overlapping genes found in GRanges object")
#    }
#    
#
#
#
#    # If there are missing assignments, print some of them
#    if (length(missing_gene_promoters) > 0) {
#        print(head(missing_gene_promoters, 10))  # Print first 10 problematic promoters
#    }
#
#    filename <- file.path(output_results_directory, glue("TC_GR_{sample}.csv"))
#    df <- as.data.frame(sample_tagclusters)
#    write.csv(df, filename, row.names = FALSE)
#    cat("Saved:", filename, "\n")
#    exit()
#}
## Save Consensus TagClusters to a csv - first need to do cumulativeCTSSdistribution, quantilePositions before calling consensusClustersGR.
##distclu_myCAGEset <- cumulativeCTSSdistribution(distclu_myCAGEset, clusters = "tagClusters", useMulticore = T)
##distclu_myCAGEset <- quantilePositions(distclu_myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
##distclu_myCAGEset <- aggregateTagClusters(distclu_myCAGEset, qLow = 0.1, qUp = 0.9, maxDist = 100)
##distclu_myCAGEset <- annotateConsensusClusters(distclu_myCAGEset, granges, upstream = 500, downstream = 500)
##
##consensus_tagclusters_GR <- consensusClustersGR(distclu_myCAGEset)
##filename <- file.path(output_results_directory, "TC_GR_Consensus.csv")
##df <- as.data.frame(consensus_tagclusters_GR)
##write.csv(df, filename, row.names = FALSE)
##cat("Saved:", filename, "\n")
##
##plot_file_path <- file.path('no_merge_reps', "TC_genomic_locations_500_500.pdf")
##plot <- plotAnnot(
##    TC,
##    scope = "counts",
##    title = "TC Genomic Location")
##ggsave(plot_file_path, plot, width = 8, height = 6, dpi = 900)
##message(glue("Saved the annotation plot to {plot_file_path}")) 
## Consensus clusters generation
##distclu_myCAGEset <- cumulativeCTSSdistribution(distclu_myCAGEset, clusters = "tagClusters", useMulticore = T)
##distclu_myCAGEset <- quantilePositions(distclu_myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
##
##distclu_myCAGEset <- aggregateTagClusters(distclu_myCAGEset, qLow = 0.1, qUp = 0.9, maxDist = 100)
##distclu_myCAGEset$outOfClusters / distclu_myCAGEset$librarySizes *100
#
##    SRO1     SRO2     SRO3     SLE1     SLE2     SLE3     SSP1     SSP2 
##41.98408 43.16077 42.49496 33.24254 33.93881 31.54148 47.04752 46.42048 
##    SSP3     SIS1     SIS2     SIS3      CR1      CR2      CR3      CR4 
##45.69298 43.08253 45.40137 46.54275 38.58052 29.23683 39.36710 38.54684 
##     CR5      CL1      CL2      CL3      CL4      CL5 
##39.78516 31.03490 38.84250 28.69027 29.71804 28.21941 
#
##distclu_myCAGEset <- annotateConsensusClusters(distclu_myCAGEset, granges, upstream = 500, downstream = 50)
#
#
#
#
##print('Export consensus')
##tagclusters <- tagClustersGR(distclu_myCAGEset)
##track_file_path <- file.path(output_results_directory, glue("Consensus_TC.bed"))
##trk <- exportToTrack(tagclusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
##export.bed(trk, con = track_file_path)
##
### Save all TagClusters to a bed file 
##for (sample in sampleLabels) {
##    # Export to track
##    message(glue("{Sys.Date()} - Exporting {sample} to track ..."))
##    sample_tagclusters <- tagClustersGR(distclu_myCAGEset, sample = sample)
##    track_file_path <- file.path(output_results_directory, glue("{sample}_TC.bed"))
##    trk <- exportToTrack(sample_tagclusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
##    export.bed(trk, con = track_file_path)
##    message(glue("{Sys.Date()} - Exported {sample} to track successfully."))
##}
###
### Extract Tag Clusters for SRO1, SRO2, and SRO3
##sro1_clusters <- tagClustersGR(distclu_myCAGEset, sample = "SRO1")
##sro2_clusters <- tagClustersGR(distclu_myCAGEset, sample = "SRO2")
##sro3_clusters <- tagClustersGR(distclu_myCAGEset, sample = "SRO3")
##merged_clusters <- reduce(c(sro1_clusters, sro2_clusters, sro3_clusters))
##merged_track_file_path <- file.path(output_results_directory, "SRO_merged_TC.bed")
##trk <- exportToTrack(merged_clusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
##export.bed(trk, con = merged_track_file_path)
##message(glue("{Sys.Date()} - Exported merged track to {merged_track_file_path} successfully."))
#
#enhancers <- quickEnhancers(distclu_myCAGEset) # enhancers: RangedSummarizedExperiment with 25606 rows and 22 columns
#enhancers_GR <- tagClustersGR(enhancers[["enhancers"]])
#
#print(enhancers)
#stop()
#enhancers <- annotateTagClusters(enhancers, granges, upstream = 500, downstream = 500)
#print(enhancers)
#track_file_path <- file.path(output_results_directory, glue("All_Enhancers.bed"))
#enhancers[enhancers$annotation %in% c("intergenic", "intron"), ]
#print(enhancers)
#enh_trk <- exportToTrack(enhancers, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
#rtracklayer::export.bed(enh_trk, track_file_path)
#stop()
#
## keep only non-exonic BCs as enhancer candidates -  subset(annotation == "intergenic" or "intron")
#
## Save all Enhancers to a csv and track 
#for (sample in sampleLabels) {
#    message(glue("{Sys.Date()} - "))
#    sample_enhancers <- tagClustersGR(enhancers, sample = sample)
#    filename <- file.path(output_results_directory, glue("Enhancer_GR_{sample}.csv"))
#    df <- as.data.frame(sample_enhancers)
#    write.csv(df, filename, row.names = FALSE)
#    cat("Saved:", filename, "\n")
#    track_file_path <- file.path(output_results_directory, glue("Enhancer_{sample}.bed"))
#    trk <- exportToTrack(sample_enhancers, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
#    export.bed(trk, con = track_file_path)
#    message(glue("{Sys.Date()} - Exported {sample} to csv and track successfully."))
#}
## Save Consensus Enhancers to a csv 
#consensus_enhancers <- consensusClustersGR(enhancers)
#filename <- file.path(output_results_directory, "Enhancer_GR_Consensus.csv")
#df <- as.data.frame(consensus_enhancers)
#write.csv(df, filename, row.names = FALSE)
#cat("Saved:", filename, "\n")
## Export to track
#track_file_path <- file.path(output_results_directory, "Enhancer_Consensus.bed")
#trk <- exportToTrack(consensus_enhancers, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
#export.bed(trk, con = track_file_path)
#message(glue("{Sys.Date()} - Exported consensus enhancers to csv and track successfully."))
#