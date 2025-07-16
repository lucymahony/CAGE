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
# Save the clustered object
output_filename <- glue("{output_results_directory}/{clusteringMethod}_md{maxDist}_single_{keepSingletonsAbove}_myCAGEset_ALL_merged.rds")
saveRDS(myCAGEset, output_filename)
message(glue("{Sys.Date()} - Clustering complete, saved to {output_filename}."))

# Export the clusters to BED format # There is only one sample, ALL
for (sample in sampleLabels(myCAGEset)) {
    message(glue("{Sys.Date()} - Exporting {sample} to track ..."))
    sample_tagclusters <- tagClustersGR(myCAGEset, sample = sample)
    track_file_path <- file.path(output_results_directory, glue("{sample}_TC.bed"))
    trk <- exportToTrack(sample_tagclusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = FALSE)
    export.bed(trk, con = track_file_path)
    message(glue("{Sys.Date()} - Exported {sample} to track successfully."))
}

# Now begining the custom annotation process
# Read in the clustered object 
#myCAGEset <- readRDS(glue("{output_results_directory}/{clusteringMethod}_md{maxDist}_single_{keepSingletonsAbove}_myCAGEset_ALL_merged.rds"))

# GFF to GRanges
convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
  print('Generating GRanges object from GFF...')
  data_frame <- read.delim(gff_file, header = FALSE, skip = 1, stringsAsFactors = FALSE,
                           colClasses = c(rep("character", 8), "character"))
  colnames(data_frame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
  data_frame$feature[data_frame$feature == "mRNA"] <- "transcript"
  data_frame <- data_frame %>%
    mutate(
      gene_name = case_when(
        grepl("ID=([^;]+)", attributes) ~ sub(".*ID=([^;]+).*", "\\1", attributes),
        grepl("gene_id=([^;]+)", attributes) ~ sub(".*gene_id=([^;]+).*", "\\1", attributes),
        grepl("Parent=([^;]+)", attributes) ~ sub(".*Parent=([^;]+).*", "\\1", attributes),
        TRUE ~ NA_character_
      )
    )
  GRanges(
    seqnames = Rle(data_frame$seqname),
    ranges = IRanges(start = as.numeric(data_frame$start), end = as.numeric(data_frame$end)),
    strand = Rle(data_frame$strand),
    type = Rle(data_frame$feature),
    transcript_type = Rle(rep(transcript_type, nrow(data_frame))),
    gene_name = Rle(data_frame$gene_name)
  )
}


granges <- convert_gff_to_granges("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")
granges_transcripts <- granges[granges$type == "transcript"]  
# Annotate CTSS and TCs
myCAGEset <- annotateCTSS(myCAGEset, granges, upstream=500, downstream=500)
myCAGEset <- annotateTagClusters(myCAGEset, granges, upstream=500, downstream=500)


# Custom nearest-gene promoter function
nearest_gene_promoter <- function(seqnames, start, end, strand, dominant_ctss_pos, granges, distance = 500) {
  tc_gr <- GRanges(seqnames = Rle(seqnames),
                   ranges = IRanges(start = as.numeric(start), end = as.numeric(end)),
                   strand = Rle(strand))
  nearest_hit <- nearest(tc_gr, granges)
  if (is.na(nearest_hit)) return("")
  nearest_gene_name <- as.character(mcols(granges)[nearest_hit, "gene_name"])
  nearest_gene_strand <- as.character(strand(granges)[nearest_hit])
  if (strand != nearest_gene_strand) return("strand_mismatch")
  if (abs(dominant_ctss_pos - start(granges)[nearest_hit]) <= distance) return(nearest_gene_name)
  return("")
}

# Process a single sample
print(glue("The sample labels are {sampleLabels(myCAGEset)}"))
for (sample in sampleLabels(myCAGEset)) {
  sample_tagcluster_gr <- tagClustersGR(myCAGEset, sample = sample)

  print(glue("Processing ...{sample}"))
  df <- as.data.frame(sample_tagcluster_gr)
  
  # Pre-filter for promoter TCs
  promoter_rows <- df$annotation == "promoter"
  promoter_df <- df[promoter_rows, ]
   
  print(glue("Promoter TCs: {nrow(promoter_df)}"))

  # Assign nearest gene only to promoter TCs
  promoter_df$nearest_gene <- Map(function(seqnames, start, end, strand, dominant_ctss_pos) {
    nearest_gene_promoter(
      seqnames = as.character(seqnames),
      start = as.numeric(start),
      end = as.numeric(end),
      strand = as.character(strand),
      dominant_ctss_pos = as.numeric(dominant_ctss_pos),
      granges = granges_transcripts
    )
  }, promoter_df$seqnames, promoter_df$start, promoter_df$end, promoter_df$strand, promoter_df$dominant_ctss.pos)

  # Post-processing
  number_tc_opposite_strand <- sum(promoter_df$nearest_gene == "strand_mismatch")
  print(glue("Opposite strand promoter TCs: {number_tc_opposite_strand}"))
  promoter_df$nearest_gene[promoter_df$nearest_gene == "strand_mismatch"] <- ""

  # Merge back
  df$nearest_gene <- ""
  df$nearest_gene[promoter_rows] <- promoter_df$nearest_gene

  # Fill in gene column where it's empty and annotation is promoter
  original_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
  df$gene <- ifelse(df$gene == "" & df$annotation == "promoter", df$nearest_gene, df$gene)
  new_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
  newly_annotated_genes <- original_empty_genes - new_empty_genes

  print(glue("Added annotations to {newly_annotated_genes} promoter TCs."))

  # Final check
  if (original_empty_genes - newly_annotated_genes != new_empty_genes) {
    stop("Mismatch in expected annotation counts.")
  }
  list_cols <- sapply(df, is.list)
  print(glue("List columns: {sum(list_cols)}"))
  if (any(list_cols)) {
    df[list_cols] <- lapply(df[list_cols], function(col) {
      sapply(col, function(x) {
        if (is.null(x)) return(NA)
        if (length(x) == 1) return(as.character(x))
        return(paste(as.character(x), collapse = ";"))
      })
    })
  }
  # Save result myCAGEset <- readRDS("{output_results_directory}/{clusteringMethod}_md{maxDist}_single_{keepSingletonsAbove}_myCAGEset_ALL_merged.rds") 
  filename <- file.path(output_results_directory, glue("{clusteringMethod}_md{maxDist}_single_{keepSingletonsAbove}_myCAGEset_ALL_merged_my_annot.csv"))
  write.csv(df, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}
