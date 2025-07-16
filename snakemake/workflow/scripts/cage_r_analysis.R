 ## The script performs the following steps 
## 1. read data in 
## 2. quality control 
## 3. merge samples 
## 4. normalise
## 5. cluster 
## 6. location_tc_s
## 6. promoter width 
## 7. consensus across samples 
## 8. shifting promoter

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
  #This function parses in the paramaters from the command line. 
  #Returns the 
  #  file pattern
  #  input files directory
  #  input files
  #  samples tsv file path

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

read_in_cage_data <- function(input_files, file_pattern, type) {
 # Reads in the ctss files and print in the library sizes 
  message(glue('{Sys.Date()} - Reading in CAGE data ....'))
  genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
  myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                       inputFiles = input_files, 
                       inputFilesType = type,
                       sampleLabels = sub(file_pattern, "", basename(input_files)))
  myCAGEset <-getCTSS(myCAGEset) 
  message(glue("{Sys.Date()} - Successfully read in CAGE files."))
  library_size <- librarySizes(myCAGEset)
  print(glue('CAGE object library size {library_size}'))
  return(myCAGEset)
}

correlation_matrix <- function(myCAGEset, output_results_directory, tagCountThreshold = 5, applyThresholdBoth = FALSE, method = "pearson") {
  #' Generate Correlation Matrix and Scatterplots
  #'
  #' This function creates a scatterplot matrix of pairwise CAGE signals and computes
  #' the correlation matrix for the samples in a given CAGEset object. The plot is 
  #' saved as a PDF file, and the correlation matrix is written to a CSV file.

  print(glue('{Sys.Date()} - Calculating correlation matrix ....'))
  if (!dir.exists(output_results_directory)) {
    stop(glue("ERROR: Output directory {output_results_directory} does not exist."))
  }
  tryCatch({
    # Generate the correlation matrix
    correlation_matrix <- plotCorrelation(
      myCAGEset, what = "CTSS", values = "raw", samples = "all",
      method = method, tagCountThreshold = tagCountThreshold, 
      applyThresholdBoth = applyThresholdBoth, plotSize = 900
    )
    # Saves the file to Rplots.pdf automatically, so we need to rename it to correlation_matrix.pdf
    old_file <- file.path(output_results_directory, "Rplots.pdf")
    new_file <- file.path(output_results_directory, "correlation_matrix.pdf")
    matrix_file <- file.path(output_results_directory, "correlation_matrix.csv")
    
    if (file.exists(old_file)) {
      if (!file.rename(from = old_file, to = new_file)) {
        stop("ERROR: Failed to rename the correlation matrix PDF file.")
      } else {
        message(glue("PDF renamed successfully to {new_file}."))
      }
    }
    write.csv(correlation_matrix, file = matrix_file, row.names = TRUE)
    print(glue('{Sys.Date()} - Correlation matrix saved to {matrix_file}, Note that the correlation matrix might appear blank before it loads as this is a large file.'))
  }, error = function(e) {
    message(glue("ERROR: Failed to calculate correlation matrix: {e$message}"))
    stop(e)
  })
}

merge_CAGE_replicates <- function(myCAGEset, samples_tsv_file_path) {
  message(glue("{Sys.Date()} - Generating labels from samples using {samples_tsv_file_path} ..."))
  samples <- tryCatch({
    # Read the samples TSV file - comment.char = '#' is important 
    read.delim(samples_tsv_file_path, header = TRUE, sep = "", stringsAsFactors = FALSE, comment.char = '#')
    }, error = function(e) {
    stop(glue("Failed to read samples TSV file: {e$message}"))
  })
  if (!all(c("sample_name", "tissue", "genome") %in% colnames(samples))) {
    stop("ERROR: The samples TSV file must contain 'sample_name', 'tissue', and 'genome' columns.")
  }
  
  samples <- samples %>%
    arrange(sample_name)
  samples <- samples %>%
    mutate(mergedLabel = paste(tissue, genome, sep = "_"))
  samples <- samples %>%
    mutate(groupNumber = match(mergedLabel, unique(mergedLabel)))
  
  mergeIndex <- samples$groupNumber
  mergedSampleLabels <- unique(samples$mergedLabel)

  myCAGEset <- tryCatch({
    mergeSamples(
      myCAGEset,
      mergeIndex = mergeIndex, 
      mergedSampleLabels = mergedSampleLabels)},
      error = function(e) {
    stop(glue("ERROR: Failed to merge samples: {e$message}"))
  })
  message(glue("{Sys.Date()} - Merging replicates completed successfully."))
  return(myCAGEset)
}

reverse_cumulatives <- function(myCAGEset, output_results_directory, method = "powerLaw", fitInRange = c(5, 1000)) {
  # Generate the reverse cumulative plot and use it to normalize the CAGEset
  message(glue("{Sys.Date()} - Generating reverse cumulative plot ..."))
  reverse_cumulative_plot <- plotReverseCumulatives(
    object = myCAGEset,
    fitInRange = fitInRange
  )
  reverse_cumulative_plot_file_path <-  file.path(output_results_directory, "Reverse_cummulatives_plot.pdf")
  ggsave(reverse_cumulative_plot_file_path, reverse_cumulative_plot, width = 8, height = 6, dpi=900)
  message(glue("{Sys.Date()} - Reverse cumulative plot saved to {reverse_cumulative_plot_file_path}"))
}



normalise_and_cluster_samples <- function(myCAGEset, output_results_directory, alpha, t, lower_range, upper_range, maxDist = 50, keepSingletonsAbove = 5) {
  # STEP 1: Normalise the tag counts
  message(glue("{Sys.Date()} - Normalising samples with an alpha of {alpha}, T of {t}, lower_range of {lower_range} and upper_range of {upper_range} ..."))
  myCAGEset <- normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(lower_range, upper_range), alpha = alpha, T = t)
  message(glue("{Sys.Date()} - Normalising finished."))
  
  # STEP 2: Cluster the samples
  message('Different clustering methods can be used, see word doc - \'CAGE_clustering_methods_lit_review\'')
  message(glue("{Sys.Date()} - Clustering samples with a maximum distance of {maxDist} and keeping singletons above {keepSingletonsAbove} ..."))
  myCAGEset <- distclu(myCAGEset, maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)
  message(glue("{Sys.Date()} - Clustering finished."))
  
  # STEP 3: Save the clustered samples as an RDS object
  saveRDS(myCAGEset, file = file.path(output_results_directory, "CAGE_object_clustered_samples.rds"))
  message(glue("{Sys.Date()} - Saved myCAGEset object to {output_results_directory}/CAGE_object_clustered_samples.rds with saveRDS."))

  # Read in CAGE from RDS
  message(glue("{Sys.Date()} - Reading in CAGE object from RDS ..."))
  myCAGEset <- readRDS(file.path(output_results_directory, "CAGE_object_clustered_samples.rds"))
  
  # STEP 4: Export clustered samples to track and save genomic ranges as a text file
  sample_names <- sampleLabels(myCAGEset)
  for (sample in sample_names) {
    # Export to track
    message(glue("{Sys.Date()} - Exporting {sample} to track ..."))
    sample_tagclusters <- tagClustersGR(myCAGEset, sample = sample)
    track_file_path <- file.path(output_results_directory, glue("{sample}_TC.bed"))
    trk <- exportToTrack(sample_tagclusters, what = "tagClusters", colorByExpressionProfile = FALSE, oneTrack = TRUE)
    export.bed(trk, con = track_file_path)
    message(glue("{Sys.Date()} - Exported {sample} to track successfully."))

    # Save tag cluster genomic ranges to a text file
    message(glue("{Sys.Date()} - Saving genomic ranges for {sample} to text file ..."))
    text_file_path <- file.path(output_results_directory, glue("{sample}_TC_width.txt"))
    write.table(sample_tagclusters, text_file_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    message(glue("{Sys.Date()} - Saved genomic ranges for {sample} to {text_file_path}."))
  }

  # STEP 5: Return the modified CAGEset object
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


plot_location_ctss <- function(myCAGEset, granges, output_results_directory) {
  print(glue("{Sys.Date()} - Plotting location of CTSS..."))
  print(glue("Contents of granges:"))
  print(granges)
  myCAGEset <- annotateCTSS(myCAGEset, granges, upstream=500, downstream=500)
  print(glue("Annotated CTSS data frame:"))
  print(colData(myCAGEset)[,c("librarySizes", "promoter", "exon", "intron", "unknown")])
  # Reorder
  custom_order <- c("SPI_chinese_spring", "SE_chinese_spring", "RO_chinese_spring",  "EM_chinese_spring",
                   "SP_fielder", "LE_fielder", "RO_fielder", "IS_fielder", 
                   "RO_cadenza", "LE_cadenza")
  col_data <- as.data.frame(colData(myCAGEset))
  col_data$group <- rownames(col_data) 
  col_data$group <- factor(col_data$group, levels = custom_order)
  message(glue("Col data 1 = {col_data}"))
  colData(myCAGEset)$group <- col_data$group  # Add the ordered group column to the CAGEexp object
  message(glue("Col data 2 = {col_data}"))
  # Plot the CTSS genomic location
  plot_file_path <- file.path(output_results_directory, "CTSS_genomic_locations_500_up_500_down.pdf")
  plot <- plotAnnot(
    myCAGEset,
    scope = "counts",
    title = "CTSS Genomic Location",
    group = "group"  # Use the custom-ordered group for the y-axis
  )
  ggsave(plot_file_path, plot, width = 8, height = 6, dpi = 900)
  message(glue("Saved the annotation plot to {plot_file_path}")) 
}


plot_iqr_tag_clusters <- function(myCAGEset, output_results_directory) {
  # This function plots the interquartile width for tag clusters - so you can see the distribution of narrow and broad promoters
  message(glue("{Sys.Date()} - cumulative CTSS distribution"))
  myCAGEset <- cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters", useMulticore = TRUE)
  message(glue("{Sys.Date()} - Quantile positions for tag clusters..."))
  myCAGEset <- quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  message(glue("{Sys.Date()} - Plotting interquartile width for tag clusters..."))
  
  iqr_plot <- plotInterquantileWidth(
    myCAGEset, 
    clusters = "tagClusters", 
    tpmThreshold = 3, 
    qLow = 0.1, 
    qUp = 0.9
  )
  plot_file <- file.path(output_results_directory, "tag_cluster_interquantile_width.pdf")
  ggsave(plot_file, iqr_plot, width = 8, height = 6, dpi = 900)
  message(glue("Interquartile width plot saved to {plot_file}"))
}


enhancers <- function(myCAGEset, output_results_directory) {
  # This function identifies enhancers in the CAGEset object
  # Enhancers are defined as regions with high tag counts that are not promoters
  sample_names <- sampleLabels(myCAGEset)
  print(glue("{Sys.Date()} - Sample names: {sample_names}"))
  for (sample in sample_names) {
    message(glue("{Sys.Date()} - Identifying enhancers for {sample}..."))
    enhancers <- quickEnhancers(myCAGEset)
    print(enhancers[["enhancers"]])
    enhancers <- enhancers[, sample]
    sample_tagclusters <- tagClustersGR(enhancers)
    trk <- exportToTrack(sample_tagclusters, what ="tagClusters",colorByExpressionProfile = FALSE,oneTrack = TRUE)
    track_file_path <- file.path(output_results_directory, glue("{sample}_enhancers.bed"))
    export.bed(trk, con=track_file_path)
    message(glue("{Sys.Date()} - Enhancers saved to {track_file_path}"))
    rds_path <- file.path(output_results_directory, glue("{sample}_enhancers.rds"))
    saveRDS(enhancers, file = rds_path)
    message(glue("{Sys.Date()} - Saved enhancers to {rds_path}"))
    }
  }
  

expression_profiling <- function(myCAGEset, output_results_directory, tpmThreshold = 10) {
  message(glue("{Sys.Date()} - Performing expression profiling..."))
  myCAGEset <- getExpressionProfiles(myCAGEset, what = 'CTSS', tpmThreshold = tpmThreshold, nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)
  expression_plot <- plotExpressionProfiles(myCAGEset, what = "CTSS")
  expression_plot_file_path <- file.path(output_results_directory, "Expression_profiles.pdf")
  ggsave(expression_plot_file_path, expression_plot, width = 8, height = 6, dpi=900)
  print('Finished expression profiling')
}

shifting_promoters <- function(myCAGEset, output_results_directory) {
  " Shifting score is a measure of differential usage of TSSs within consensus cluster between two samples"
  print('Consensus clusters must be calculated before shifting promoters')
  shifting.promoters <- getShiftingPromoters(myCAGEset, 
    groupX = "CL1", groupY = "CL2",
        tpmThreshold = 5, scoreThreshold = 0.6,
        fdrThreshold = 0.01)
  head(shifting.promoters)
  print('Finished promoter shifting')
}


annotated_ctss <- function(myCAGEset, granges, output_results_directory, csv_file_name, upstream = 500, downstream = 500) {
  # Annotate the CTSS data
  annotated_object <- annotateCTSS(myCAGEset, granges, upstream = upstream, downstream = downstream)
  print(glue("Annotated CTSS data frame:"))
  print(colData(annotated_object)[,c("librarySizes", "promoter", "exon", "intron", "unknown")])
  annotated_ctss <- CTSScoordinatesGR(annotated_object)
  # Convert to GRanges with Metadata
  annotated_gr <- GRanges(
    seqnames = seqnames(annotated_ctss),
    ranges = ranges(annotated_ctss),
    strand = strand(annotated_ctss),
    annotation = mcols(annotated_ctss)$annotation, # Promoter, Exon, etc.
    genes = mcols(annotated_ctss)$genes           # Associated gene names
  )
  file_name <- glue("{output_results_directory}/up_{upstream}_down_{downstream}_{csv_file_name}")
  write.table(annotated_gr, file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

annotating_genes <- function(list_input_files, gff_file, upstream, downstream, output_results_directory, type) {
  # For prior generated TagClusters, they are read in using existing function read_in_cage_data pretending the tag clusters are 'ctss'
  # And saves the annotated LC genes to a csv file for each tissue/genome TagCluster file
  granges <- convert_gff_to_granges(gff_file, transcript_type = "protein_coding")
  for (file in list_input_files) {
    message(glue("{Sys.Date()} - Annotating {type} genes for file {file}"))
    file_pattern <-"_at_least_2_overlaps_max_score.bed"
    myCAGEset <- read_in_cage_data(file, file_pattern, "bed")
    # Name the csv file after the input files
    csv_file_name <- sub(file_pattern, "", basename(file))
    csv_file_name <- paste0(csv_file_name, "_annotated_", type, "_genes_ctss.csv")
    annotated_ctss(myCAGEset, granges, output_results_directory, csv_file_name, upstream, downstream)
  }
}

count_genes_with_tag_clusters <- function(list_input_files, output_results_directory, upstream, downstream, type) {
  # Uses the annotated LC genes csv files to count the number of LC genes that are supported by TC's.
  gene_tissue_map <- list()
  for (file in list_input_files) {
    tissue_name <- sub("_at_least_2_overlaps_max_score.bed", "", basename(file))
    message(glue("{Sys.Date()} - Counting {type} genes with tag clusters for file {file}"))
    csv_file_name <- paste0("/up_", upstream, "_down_", downstream, "_", tissue_name, "_annotated_", type, "_genes_ctss.csv")
    annotated_ctss_data <- read.table(file = paste0(output_results_directory, csv_file_name), header = TRUE, sep = "\t")
    # Subset the data frame to only include the promoter annotation
    annotated_ctss_data_promoter <- annotated_ctss_data[annotated_ctss_data$annotation == "promoter",]
    gene_list <- unique(annotated_ctss_data_promoter$genes)
    print(glue("Number of unique genes = {length(gene_list)}"))
    for (gene in gene_list) {
      if (!gene %in% names(gene_tissue_map)) {
        gene_tissue_map[[gene]] <- list()
      }
      gene_tissue_map[[gene]] <- unique(c(gene_tissue_map[[gene]], tissue_name))
    }
  }
  gene_names <- names(gene_tissue_map)
  tissues_present_in <- sapply(gene_tissue_map, function(tissues) paste(tissues, collapse = ", "))
  gene_tissue_df <- data.frame(Gene_name = gene_names, Tissues_present_in = tissues_present_in, stringsAsFactors = FALSE)
  gene_tissue_df <- gene_tissue_df[!is.na(gene_tissue_df$Gene_name),]
  gene_tissue_df$Gene_name <- gsub("^NA;", "", gene_tissue_df$Gene_name)
  gene_tissue_df <- gene_tissue_df[gene_tissue_df$Gene_name != "",]
  write.csv(gene_tissue_df, paste0(output_results_directory, "gene_tissue_map_lc.csv"), row.names = FALSE)
  # Print key numbers 
  genes_in_all_tissues <- sum(sapply(gene_tissue_map, length) == length(list_input_files))
  genes_in_one_tissue <- sum(sapply(gene_tissue_map, length) == 1)
  total_genes <- length(gene_tissue_map)
  print(glue("Number of genes present in all tissues: {genes_in_all_tissues}"))
  print(glue("Number of genes present in only one tissue: {genes_in_one_tissue}"))
  print(glue("Total number of genes that have support: {total_genes}"))  
  return(gene_tissue_df)
}




main <- function() {
  # STEP 1: Define file paths and convert GFF to GRanges
  gff_file_path_hc <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
  #hc_granges <- convert_gff_to_granges(gff_file_path_hc) 
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
  
  # STEP 2: Parse command-line arguments, read in cage data, and set ouput directory
  args <- parse_arguments()  
  #myCAGEset <- read_in_cage_data(args$input_files, args$file_pattern, "ctss") 
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
  #
  ## STEP 3: Generate correlation matrix 
  #correlation_matrix(myCAGEset, output_results_directory, tagCountThreshold = 5, applyThresholdBoth = FALSE, method = "pearson") 

  ## STEP 4: - Merging replicates and reverse cumulative plot generation (commented out for now)
  #myCAGEset <- merge_CAGE_replicates(myCAGEset, args$samples_tsv_file_path)

  ## STEP 5: Reverse cumulative plot generation - to get parameters 
  #reverse_cumulatives(myCAGEset, output_results_directory, method = "powerLaw", fitInRange = c(5, 1000))
  #saveRDS(myCAGEset, file = "myCAGEset_before_normalisation_and_clustering.rds")
  #alpha = 1.06
  #t = 1e+07
  #lower_range = 1
  #upper_range = 100000

  ## STEP 6: Normalise and cluster samples
  #myCAGEset <- normalise_and_cluster_samples(myCAGEset, output_results_directory, alpha, t, lower_range, upper_range, maxDist = 20, keepSingletonsAbove = 5)

  ## load the CAGEset object after normalisation and clustering
  myCAGEset <- readRDS("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/myCAGEset_with_enhancers.rds")
  ## STEP 7: Enhancer identification
  #enhancers(myCAGEset, output_results_directory)

  ## STEP 9: Plot data 
  #plot_location_ctss(myCAGEset, hc_granges, output_results_directory)
  #plot_iqr_tag_clusters(myCAGEset, output_results_directory)

  ## STEP 12: Expression profiling
  #expression_profiling(myCAGEset, output_results_directory, tpmThreshold = 5)

  # Step 13: Annotate HC and LC genes
  ## List the input files as all the files in the output directory that end in   at_least_2_overlaps_max_score.bed
  gff_file_paths <- list(
  hc = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3",
  lc = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_LC.gff3")
  list_input_files <- list.files("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data", pattern = "_at_least_2_overlaps_max_score.bed", full.names = TRUE)
  for (type in names(gff_file_paths)) {
    gff_file_path <- gff_file_paths[[type]]
    annotating_genes(list_input_files, gff_file_path, 500, 500, output_results_directory, type)
    gene_tissue_df <- count_genes_with_tag_clusters(list_input_files, output_results_directory, 500, 500, type)
  }
  ## Additional steps
  #Calculate the shifting promoters and then 
  #shifting_promoters(myCAGEset, output_results_directory)

#
}


main()

