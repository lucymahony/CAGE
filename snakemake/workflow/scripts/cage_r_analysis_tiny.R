# This is a tiny copy of CAGE r so I can mess around with aggregate tag clusters and consensus clusters. 
suppressPackageStartupMessages({
  library(CAGEr)
  library(BSgenome) # If there are any issues with BSgenome check that all the packages have been updated and compiled
  library(BSgenome.Taestivum.ChineseSpring)
  library(rtracklayer)
  library(GenomicRanges)
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(corrplot)
  library(tibble)
  library(base)
  library(glue)
})


parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # Validate thre arguments are provided
  if (length(args) < 3) {
    stop(" ERROR: Please provide three commandline arguments: <file_pattern> <input_files_directory> <samples_tsv_file_path>")
  }
  file_pattern <- args[1]
  input_files_directory <- args[2]
  input_files <- list.files(path = input_files_directory, pattern = file_pattern, full.names = TRUE)
  # Check if all listed files exist
  file_existence <- file.exists(input_files)
  if (!all(file_existence)) {
    missing_files <- input_files[!file_existence]
    stop(
      "The following files do not exist:\n", 
      paste(missing_files, collapse = "\n")
    )
  } else {
    message("All files exist.")
    message(file_existence)
  }
  samples_tsv_file_path <- args[3]
  message("WARNING: The successful running of this script requires the config files to be correctly formatted.")
  message("Check the samples TSV file contains the samples you wish to process ")

  # Return the parsed arguments and input files as named lists
  list(
    file_pattern = file_pattern,
    input_files_directory = input_files_directory,
    input_files = input_files,
    samples_tsv_file_path = samples_tsv_file_path
  )
}

read_in_cage_data <- function(input_files, sampleLabels, type) {
 # Reads in the ctss files and print in the library sizes 
  message(glue('{Sys.Date()} - Reading in CAGE data ....'))
  genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
  myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                       inputFiles = input_files, 
                       inputFilesType = type,
                       sampleLabels = sampleLabels)
  myCAGEset <-getCTSS(myCAGEset) 
  message(glue("{Sys.Date()} - Successfully read in CAGE files."))
  library_size <- librarySizes(myCAGEset)
  print(glue('CAGE object library size {library_size}'))
  return(myCAGEset)
}

convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
  print('Generating own granges object from GFF file .... ')
  #print('Check that generated granges object has the same structure as the exampleZv9_annot object')
  #print('Head of exampleZv9_annot = ')
  #print(head(exampleZv9_annot))

  data_frame <- read.delim(gff_file, header=F, skip = 1)
  data_frame <- data_frame %>%
    mutate(extracted_id = ifelse(
    grepl('=([^\\.]+)\\.', V9),
    sub('.*=([^\\.]+)\\..*', '\\1', V9),
    NA)
  )
  # Replace "mRNA" with "transcript" - Thiis is important for the plot annot to work!
  data_frame$V3[data_frame$V3 == "mRNA"] <- "transcript"
  seqnames <- Rle(data_frame$V1) # e.g. Chr1A 
  targetRanges <- IRanges(data_frame$V4, data_frame$V5) 
  strands <- Rle(data_frame$V7)
  gene_name <- Rle(data_frame$extracted_id)
  type <- Rle(data_frame$V3)
  
  transcript_type <- Rle(rep(transcript_type, length(seqnames)))
  gr <- GRanges(seqnames = seqnames, ranges = targetRanges, strand = strands, type = type, transcript_type = transcript_type, gene_name = gene_name)
  return(gr)
}



read_to_aggregate_TC <- function(input_files, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path) {
  gff_file_path_hc <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
  
  print(input_files)
  ce <- read_in_cage_data(input_files, sampleLabels, "ctss") 
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
  print('merging....')
  ce <- mergeSamples(ce, mergeIndex = mergeIndex, mergedSampleLabels = mergedSampleLabels)
  print('normalizing....')
  ce <- normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
  print('filtering....')
  ce <- filterLowExpCTSS(ce, thresholdIsTpm = TRUE, nrPassThreshold = 1, threshold = 1)
  print('clustering....')
  # distclu
  ce <- distclu(ce, maxDist = 20, keepSingletonsAbove = 5)
  print('quantile positions')
  ce <- cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
  ce <- quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  print('aggregateTagClusters....')
  ce <- aggregateTagClusters(ce, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
  print(ce$outOfClusters)
  print(ce$librarySize)
  print(ce$outOfClusters / ce$librarySizes *100)
  consensusClustersGR(ce)

  # Annotate
  gff_file_path_hc <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
  granges <- convert_gff_to_granges(gff_file_path_hc) 
  print(glue("Contents of granges:"))
  print(granges)
  myCAGEset <- annotateCTSS(ce, granges, upstream=500, downstream=500)

  ce <- cumulativeCTSSdistribution(ce, clusters = "consensusClusters", useMulticore = TRUE)
  ce <- quantilePositions(ce, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)
  CCGR <- consensusClustersGR(ce, qLow = 0.1, qUp = 0.9)

  trk <- exportToTrack(CCGR, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = FALSE)
  export.bed(trk, con = track_file_path)
  message(glue("{Sys.Date()} - Exported to track successfully."))
}


main <- function() {
  ## Cadenza
  input_files_cadenza <- c('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL1.CS.se.star.unique.sorted.ctss.n.bed',
   '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL2.CS.se.star.unique.sorted.ctss.n.bed',
   '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL3.CS.se.star.unique.sorted.ctss.n.bed',
   '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL4.CS.se.star.unique.sorted.ctss.n.bed',
   '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL5.CS.se.star.unique.sorted.ctss.n.bed',
    '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR1.CS.se.star.unique.sorted.ctss.n.bed',
    '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR2.CS.se.star.unique.sorted.ctss.n.bed',
    '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR3.CS.se.star.unique.sorted.ctss.n.bed',
    '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed',
    '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR5.CS.se.star.unique.sorted.ctss.n.bed')
  sampleLabels <- c('CL1', 'CL2', 'CL3', 'CL4', 'CL5', 'CR1', 'CR2', 'CR3', 'CR4', 'CR5')
  mergeIndex <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
  mergedSampleLabels <- c("CL", "CR")
  track_file_path <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/aggregate_tag_clusters_cadenza_false.bed"
  read_to_aggregate_TC(input_files_cadenza, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path)
  ## Fielder
  #input_files_fielder <- c('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS3.CS.se.star.unique.sorted.ctss.n.bed')
  #sampleLabels <- c('SLE1', 'SLE2', 'SLE3', 'SRO1', 'SRO2', 'SRO3', 'SSP1', 'SSP2', 'SSP3', 'SIS1', 'SIS2', 'SIS3')
  #mergeIndex <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
  #mergedSampleLabels <- c("LE", "RO", "SP", "IS")
  #track_file_path <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/aggregate_tag_clusters_fielder.bed"
  #read_to_aggregate_TC(input_files_fielder, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path)
  ## Chinese Spring  
  #input_file_chinese_spring <- c(
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSRO1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSRO2.CS.se.star.unique.sorted.ctss.n.bed',
#
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE3R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE3R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE4.CS.se.star.unique.sorted.ctss.n.bed',
#
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM2R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM2R2.CS.se.star.unique.sorted.ctss.n.bed',
#
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI2R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI2R2.CS.se.star.unique.sorted.ctss.n.bed')
  #sampleLabels <-c('CSRO1', 'CSRO2', 'CSSE1R1', 'CSSE1R2', 'CSSE2', 'CSSE3R1', 'CSSE3R2', 'CSSE4', 'CSEM1R1', 'CSEM1R2', 'CSEM2R1', 'CSEM2R2', 'CSSPI1R1', 'CSSPI1R2', 'CSSPI2R1', 'CSSPI2R2')
  #mergeIndex <- c(1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
  #mergedSampleLabels <- c("RO", "SE", "EM", "SPI")
  #track_file_path <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/aggregate_tag_clusters_chinese_spring.bed"
  #read_to_aggregate_TC(input_file_chinese_spring, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path)
  ## Fielder and Cadenza
  #input_files_fielder_and_cadenza <- c('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL1.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL2.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL3.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL4.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL5.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR1.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR2.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR3.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR5.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS3.CS.se.star.unique.sorted.ctss.n.bed')
  #sampleLabels <- c('CL1', 'CL2', 'CL3', 'CL4', 'CL5', 'CR1', 'CR2', 'CR3', 'CR4', 'CR5', 'SLE1', 'SLE2', 'SLE3', 'SRO1', 'SRO2', 'SRO3', 'SSP1', 'SSP2', 'SSP3', 'SIS1', 'SIS2', 'SIS3')
  #mergeIndex <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6)
  #mergedSampleLabels <- c("CL", "CR", "LE", "RO", "SP", "IS")
  #track_file_path <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/aggregate_tag_clusters_fielder_and_cadenza.bed"
  #read_to_aggregate_TC(input_files_fielder_and_cadenza, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path)
  # Fielder cadenza and chinese spring
  #input_files_fielder_and_cadenza_and_chinese_spring <-  c('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL1.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL2.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL3.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL4.CS.se.star.unique.sorted.ctss.n.bed',
  # '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CL5.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR1.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR2.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR3.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CR5.CS.se.star.unique.sorted.ctss.n.bed',
  #  '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SLE3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SRO3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SSP3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/SIS3.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSRO1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSRO2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE3R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE3R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSE4.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM2R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSEM2R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI1R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI1R2.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI2R1.CS.se.star.unique.sorted.ctss.n.bed',
  #'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CSSPI2R2.CS.se.star.unique.sorted.ctss.n.bed')

  #sampleLabels <- c('CL1', 'CL2', 'CL3', 'CL4', 'CL5', 'CR1', 'CR2', 'CR3', 'CR4', 'CR5', 'SLE1', 'SLE2', 'SLE3', 'SRO1', 'SRO2', 'SRO3', 'SSP1', 'SSP2', 'SSP3', 'SIS1', 'SIS2', 'SIS3', 'CSRO1', 'CSRO2', 'CSSE1R1', 'CSSE1R2', 'CSSE2', 'CSSE3R1', 'CSSE3R2', 'CSSE4', 'CSEM1R1', 'CSEM1R2', 'CSEM2R1', 'CSEM2R2', 'CSSPI1R1', 'CSSPI1R2', 'CSSPI2R1', 'CSSPI2R2')
  #mergeIndex <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8,9, 9, 9, 9, 10, 10, 10, 10)
  #mergedSampleLabels <- c("CL", "CR", "LE", "RO", "SP", "IS", "RO", "SE", "EM", "SPI")
  #track_file_path <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/aggregate_tag_clusters_fielder_and_cadenza_and_chinese_spring.bed"
  #read_to_aggregate_TC(input_files_fielder_and_cadenza_and_chinese_spring, sampleLabels, mergeIndex, mergedSampleLabels, track_file_path)
}
main()

