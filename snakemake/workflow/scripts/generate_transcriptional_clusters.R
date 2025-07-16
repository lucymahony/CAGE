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
  input_files <- c("SRO1.CS.se.star.bam", "SRO2.CS.se.star.bam", "SRO3.CS.se.star.bam", "CR1.CS.se.star.bam", "CR2.CS.se.star.bam", "CR3.CS.se.star.bam", "CR4.CS.se.star.bam", "CR5.CS.se.star.bam")
  sampleLabels <- c("SRO1", "SRO2", "SRO3", "CR1", "CR2", "CR3", "CR4", "CR5")
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
# Normalize to TPM
message(glue("{Sys.Date()} - Normalising Tag Count ...."))
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
message(glue("{Sys.Date()} - Clustering...."))

if (clusteringMethod == "distclu") {
  myCAGEset <- distclu(myCAGEset, maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)
} else if (clusteringMethod == "paraclu") {
  myCAGEset <- paraclu(myCAGEset, maxLength = maxDist, keepSingletonsAbove = keepSingletonsAbove, useMulticore = TRUE, nrCores = cores)
} else {
  stop("Invalid clustering method. Must be either 'distclu' or 'paraclu'.")
}

output_filename <- glue("{output_results_directory}/{clusteringMethod}_md{maxDist}_single_{keepSingletonsAbove}_myCAGEset.rds")
saveRDS(myCAGEset, output_filename)
message(glue("{Sys.Date()} - Clustering complete, saved to {output_filename}."))


