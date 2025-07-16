## Load libraries
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
  library(BiocParallel)
})


# Get command-line arguments - Defaults are set 
args <- commandArgs(trailingOnly = TRUE)
data_directory <- ifelse(length(args) >= 1 && args[1] %in% c("merge_reps", "no_merge_reps"), args[1], "merge_reps") 
input_rds_file_name <- ifelse(length(args) >= 2, args[2], "distclu_myCAGEset.rds")
input_base_name <- tools::file_path_sans_ext(basename(input_rds_file_name))
sample_labels <-ifelse(length(args) >= 3, args[3], "SRO2,SRO3,SLE1,SLE2,SLE3,SSP1,SSP2,SSP3,SIS1,SIS2,SIS3,CR1,CR2,CR3,CR4,CR5,CL1,CL2,CL3,CL4,CL5")
sample_labels <- unlist(strsplit(sample_labels, ","))
print(glue("Data directory: {data_directory}"))
print(glue("Input RDS file name: {input_rds_file_name}"))
print(glue("Sample labels: {sample_labels}"))
print(glue("Input base name: {input_base_name}"))


# Set working dirs and load data
snakemake_intermediate_data <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
setwd(snakemake_intermediate_data)
output_results_directory <- file.path(snakemake_intermediate_data, data_directory)
distclu_myCAGEset <- readRDS(file.path(output_results_directory, input_rds_file_name)) # e.g. merge_reps/distclu_myCAGEset.rds


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
distclu_myCAGEset <- annotateCTSS(distclu_myCAGEset, granges, upstream=500, downstream=500)
distclu_myCAGEset <- annotateTagClusters(distclu_myCAGEset, granges, upstream=500, downstream=500)


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
print(glue("The sample labels are {sampleLabels(distclu_myCAGEset)}"))
for (sample in sample_labels) {
  sample_tagcluster_gr <- tagClustersGR(distclu_myCAGEset, sample = sample)

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
  # Save result
  filename <- file.path(output_results_directory, glue("TC_GR_{sample}_{input_base_name}_my_annot.csv"))
  write.csv(df, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}


