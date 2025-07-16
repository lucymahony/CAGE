# This script merges biological reps, and then clusters the merged data before creating aggregate clusters 

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

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Set defaults in case no args are passed
clusteringMethod <- ifelse(length(args) >= 1 && args[1] %in% c("distclu", "paraclu"), args[1], "distclu")
maxDist <- ifelse(length(args) >= 2 && !is.na(as.numeric(args[2])), as.numeric(args[2]), 50)
keepSingletonsAbove <- ifelse(length(args) >= 3 && !is.na(as.numeric(args[3])), as.numeric(args[3]), 5)
merge_strategy <- ifelse(length(args) >= 4 && args[4] %in% c("no_merge", "merge"), args[4], "no_merge") 
cores <- ifelse(length(args) >= 5 && !is.na(as.numeric(args[5])), as.numeric(args[5]), 1) 


message(glue("Using clustering method {clusteringMethod} with maxDist = {maxDist}, keepSingletonsAbove = {keepSingletonsAbove}, with merge strategy {merge_strategy}, and {cores} cores."))
setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data")
output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps"
myCAGEset <- readRDS('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/myCAGEset_correctsystematicG.rds')
# Merge biological replicates
message(glue("{Sys.Date()} - Merging biological replicates ...."))
#myCAGEset <- mergeSamples(myCAGEset, mergeIndex = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6), mergedSampleLabels = c("SIS", "SLE", "SRO", "SSP", "CL", "CR"))
# To speed up the process, we can use a subset of the data for testing - SIS and SLE
myCAGEset <- mergeSamples(myCAGEset, mergeIndex = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), mergedSampleLabels = c("ALL"))


# Normalize to TPM
message(glue("{Sys.Date()} - Normalising Tag Count ...."))
myCAGEset <- normalizeTagCount(myCAGEset, method =  "simpleTpm")
message(glue("{Sys.Date()} - Clustering...."))

if (clusteringMethod == "distclu") {
  myCAGEset <- distclu(myCAGEset, maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)
} else if (clusteringMethod == "paraclu") {
  myCAGEset <- paraclu(myCAGEset, maxLength = maxDist, keepSingletonsAbove = keepSingletonsAbove, useMulticore = TRUE, nrCores = cores)
} else {
  stop("Invalid clustering method. Must be either 'distclu' or 'paraclu'.")
}

convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
  print('Generating own GRanges object from GFF file .... ')
  data_frame <- read.delim(gff_file, header = FALSE, skip = 1, stringsAsFactors = FALSE, colClasses = c(rep("character", 8), "character"))
  colnames(data_frame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  data_frame$feature[data_frame$feature == "mRNA"] <- "transcript"
  if (!is.character(data_frame$attributes)) {
    data_frame$attributes <- as.character(data_frame$attributes)
  }
  data_frame <- data_frame %>%
    mutate(
      gene_name = case_when(
        grepl("ID=([^;]+)", attributes) ~ sub(".*ID=([^;]+).*", "\\1", attributes),
        grepl("gene_id=([^;]+)", attributes) ~ sub(".*gene_id=([^;]+).*", "\\1", attributes),
        grepl("Parent=([^;]+)", attributes) ~ sub(".*Parent=([^;]+).*", "\\1", attributes),
        TRUE ~ NA_character_  # Assign NA if no match is found
      )
    )
  gr <- GRanges(
    seqnames = Rle(data_frame$seqname),
    ranges = IRanges(start = as.numeric(data_frame$start), end = as.numeric(data_frame$end)),
    strand = Rle(data_frame$strand),
    type = Rle(data_frame$feature),
    transcript_type = Rle(rep(transcript_type, nrow(data_frame))),
    gene_name = Rle(data_frame$gene_name)
  )
}

granges <- convert_gff_to_granges("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")

# Export the clusters - already completed. 
for (sample in sampleLabels(myCAGEset)) {
    message(glue("{Sys.Date()} - Exporting {sample} to track ..."))
    sample_tagclusters <- tagClustersGR(myCAGEset, sample = sample)
    track_file_path <- file.path(output_results_directory, glue("{sample}_TC.bed"))
    trk <- exportToTrack(sample_tagclusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = FALSE)
    export.bed(trk, con = track_file_path)
    message(glue("{Sys.Date()} - Exported {sample} to track successfully."))
}
# Aggregate tag clusters

# Annotate 
myCAGEset <- annotateTagClusters(myCAGEset, granges)
filename <- file.path(output_results_directory, glue("TC_GR_all_.csv"))
df <- as.data.frame(tagClustersGR(myCAGEset)) 
write.csv(all, filename, row.names = FALSE)
cat("Saved:", filename, "\n")
