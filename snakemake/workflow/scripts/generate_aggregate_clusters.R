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
if (merge_strategy =="no_merge"){
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps"
  input_files <- c(
    "SIS1.CS.se.star.bam", "SIS2.CS.se.star.bam", "SIS3.CS.se.star.bam", 
    "SLE1.CS.se.star.bam", "SLE2.CS.se.star.bam", "SLE3.CS.se.star.bam",
    "SRO1.CS.se.star.bam", "SRO2.CS.se.star.bam", "SRO3.CS.se.star.bam",
    "SSP1.CS.se.star.bam", "SSP2.CS.se.star.bam", "SSP3.CS.se.star.bam", 
    "CL1.CS.se.star.bam", "CL2.CS.se.star.bam", "CL3.CS.se.star.bam", "CL4.CS.se.star.bam", "CL5.CS.se.star.bam",
    "CR1.CS.se.star.bam", "CR2.CS.se.star.bam", "CR3.CS.se.star.bam", "CR4.CS.se.star.bam", "CR5.CS.se.star.bam")
  sampleLabels <- c(
    "SIS1", "SIS2", "SIS3", 
    "SLE1", "SLE2", "SLE3",
    "SRO1", "SRO2", "SRO3",
    "SSP1", "SSP2", "SSP3", 
    "CL1", "CL2", "CL3", "CL4", "CL5",
    "CR1", "CR2", "CR3", "CR4", "CR5")
} else if (merge_strategy == "merge"){
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/merge_reps"
  input_files <- c("SRO_merged.bam", "CR_merged.bam")
  sampleLabels <- c("SRO", "CR")
} else {
  stop("Invalid merge strategy. Must be either 'no_merge' or 'merge'.")
}

message(glue("Input files: {input_files}"))

genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                     inputFiles = input_files, 
                     inputFilesType = "bam",
                     sampleLabels = sampleLabels)
message(glue("{Sys.Date()} - correctSystematicG ...."))
myCAGEset <- getCTSS(myCAGEset, correctSystematicG = TRUE)
myCAGEset <- readRDS('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/myCAGEset_correctsystematicG.rds')
# Merge biological replicates
message(glue("{Sys.Date()} - Merging biological replicates ...."))
myCAGEset <- mergeSamples(myCAGEset, mergeIndex = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6), mergedSampleLabels = c("SIS", "SLE", "SRO", "SSP", "CL", "CR"))
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

myCAGEset <- cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
myCAGEset <- quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
myCAGEset <- aggregateTagClusters(myCAGEset, tpmThreshold = 0.3, maxDist = 100, qLow = NULL, qUp = NULL)


consensus_rds_filename <- glue("consensusTCs_{clusteringMethod}_md{maxDist}.rds")
saveRDS(myCAGEset, consensus_rds_filename)

# --- Export consensus clusters to BED ---
consensus_gr <- consensusClustersGR(myCAGEset)

# Sanitize score
raw_score <- as.numeric(consensus_gr$score)
raw_score[is.na(raw_score)] <- 0
scaled_score <- round(1000 * raw_score / max(raw_score, na.rm = TRUE))

bed_consensus <- GRanges(
  seqnames = seqnames(consensus_gr),
  ranges = ranges(consensus_gr),
  strand = strand(consensus_gr),
  name = names(consensus_gr),
  score = scaled_score
)

consensus_bed_file <- glue("consensus_TCs_{clusteringMethod}_md{maxDist}.bed")
export.bed(bed_consensus, consensus_bed_file)
message(glue("Exported consensus TCs to: {consensus_bed_file}"))

# --- Export dominant CTSS positions to BED ---
dominant_gr <- consensus_gr$dominant_ctss

bed_dominant <- GRanges(
  seqnames = seqnames(dominant_gr),
  ranges = IRanges(start = start(dominant_gr), end = start(dominant_gr) + 1),
  strand = strand(dominant_gr),
  name = names(consensus_gr),
  score = scaled_score
)

dominant_bed_file <- glue("consensus_dominant_TCs_{clusteringMethod}_md{maxDist}.bed")
export.bed(bed_dominant, dominant_bed_file)
message(glue("Exported dominant CTSSs to: {dominant_bed_file}"))

# --- Generate TPM and presence/absence matrices ---
samples <- sampleLabels(myCAGEset)
consensus_ids <- names(consensus_gr)

tpm_matrix <- matrix(0, nrow = length(consensus_ids), ncol = length(samples))
rownames(tpm_matrix) <- consensus_ids
colnames(tpm_matrix) <- samples

for (i in seq_along(samples)) {
  sample <- samples[i]
  sample_cc <- consensusClustersGR(myCAGEset, sample = sample, qLow = 0.1, qUp = 0.9)
  overlapping_ids <- names(sample_cc)
  tpm_matrix[overlapping_ids, i] <- sample_cc$tpm
}

# Save TPM matrix
tpm_csv <- glue("consensus_TPM_matrix_{clusteringMethod}_md{maxDist}.csv")
write.csv(tpm_matrix, tpm_csv)
message(glue("Saved TPM matrix to: {tpm_csv}"))
# Generate presence/absence matrix
presence_matrix <- ifelse(tpm_matrix > 0, 1, 0)
presence_csv <- glue("consensus_presence_matrix_{clusteringMethod}_md{maxDist}.csv")
write.csv(presence_matrix, presence_csv)
message(glue("Saved presence/absence matrix to: {presence_csv}"))