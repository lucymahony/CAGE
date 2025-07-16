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
input_files <- c("SRO1.CS.se.star.bam", "SRO2.CS.se.star.bam", "SRO3.CS.se.star.bam", "CR1.CS.se.star.bam", "CR2.CS.se.star.bam", "CR3.CS.se.star.bam", "CR4.CS.se.star.bam", "CR5.CS.se.star.bam")
sampleLabels <- c("SRO1", "SRO2", "SRO3", "CR1", "CR2", "CR3", "CR4", "CR5")
genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")

myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                     inputFiles = input_files, 
                     inputFilesType = "bam",
                     sampleLabels = sampleLabels)
message(glue("{Sys.Date()} - correctSystematicG ...."))
myCAGEset <- getCTSS(myCAGEset, correctSystematicG = TRUE, useMulticore = TRUE, nrCores = 16)
message(glue("{Sys.Date()} - Merging samples ...."))
myCAGEset <- mergeSamples(myCAGEset, mergeIndex = c(1, 1, 1, 2, 2, 2, 2, 2), 
                   mergedSampleLabels = c("SRO", "CR"))

message(glue("{Sys.Date()} - Normalising Tag Count ...."))
myCAGEset <- normalizeTagCount(myCAGEset, method =  "simpleTpm")
maxDist <- 20
keepSingletonsAbove <- 5
message(glue("{Sys.Date()} - Clustering with the parameters: maxDist = {maxDist}, keepSingletonsAbove = {keepSingletonsAbove} ...."))
myCAGEset <- distclu(myCAGEset, maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)


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

message(glue("{Sys.Date()} - Annotating ...."))
myCAGEset <- annotateCTSS(myCAGEset, granges, upstream=500, downstream=500)
annotate <- annotateTagClusters(myCAGEset, granges, upstream = 500, downstream = 500)
message(glue("{Sys.Date()} - Cumulative distribution ...."))
myCAGEset <- cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters", useMulticore = TRUE, nrCores = 16)
message(glue("{Sys.Date()} - Score shift ...."))
myCAGEset <- scoreShift(myCAGEset, groupX = "SRO", groupY = "CR",
        testKS = TRUE, useTpmKS = FALSE)
message(glue("{Sys.Date()} - Shifting promoters ...."))
shifting.promoters <- getShiftingPromoters(myCAGEset, 
    groupX = "SRO", groupY = "CR",
        tpmThreshold = 5, scoreThreshold = 0.6,
        fdrThreshold = 0.01)
head(shifting.promoters)
message(glue("{Sys.Date()} - Saving shifting promoters ...."))
output_filename <- glue("test_shifting_promoters_SRO_CR.rds")
# Shifting.promoters is a dataframe save it as a csv
# print data type of shifting.promoters

write.csv(shifting.promoters, output_filename)
message(glue("{Sys.Date()} - Shifting promoters saved to {output_filename}."))
print(class(shifting.promoters))