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
sampleLabels <- c("SRO1", "SRO2", "SRO3", "SLE1", "SLE2", "SLE3", "SSP1", "SSP2", "SSP3", "SIS1", "SIS2", "SIS3", "CR1", "CR2", "CR3", "CR4", "CR5", "CL1", "CL2", "CL3", "CL4", "CL5")
#distclu_myCAGEset <- readRDS("no_merge_reps/distclu_myCAGEset_custom_annotated_promoters.rds")
distclu_myCAGEset <- readRDS("no_merge_reps/distclu_myCAGEset.rds")
# Save TCs to csv. 

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
distclu_myCAGEset <- annotateCTSS(distclu_myCAGEset, granges, upstream=500, downstream=500)
print(distclu_myCAGEset)

annotate <- annotateTagClusters(distclu_myCAGEset, granges, upstream = 500, downstream = 500)
for (sample in sampleLabels) {
    message(glue("{Sys.Date()} - "))
    sample_tagclusters <- tagClustersGR(annotate, sample = sample)
    missing_gene_promoters <- sample_tagclusters[sample_tagclusters$annotation == "promoter" & sample_tagclusters$genes == ""]
    # Print summary
    print(glue("Total promoters without transcript ID: {length(missing_gene_promoters)}"))

    filename <- file.path(output_results_directory, glue("TC_GR_{sample}.csv"))
    df <- as.data.frame(sample_tagclusters)
    write.csv(df, filename, row.names = FALSE)
    cat("Saved:", filename, "\n")
    stop()
}