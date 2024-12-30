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

read_in_cage_data <- function(input_files, file_pattern) {

 # Reads in the ctss files and print in the library sizes 

  message(glue('{Sys.Date()} - Reading in CAGE data ....'))
  genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
  chromosomes <- seqnames(genome)
  chrom_lengths <- seqlengths(genome)
  chrom_lengths


  myCAGEset <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                       inputFiles = input_files, 
                       inputFilesType = "ctss",
                       sampleLabels = sub(file_pattern, "", basename(input_files)))

  myCAGEset <-getCTSS(myCAGEset) 
  warnings()
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
  #'
  print(glue('{Sys.Date()} - Calculating correlation matrix ....'))
  # Ensure output directory exists
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
    
    # Save the correlation matrix
    write.csv(correlation_matrix, file = matrix_file, row.names = TRUE)
    print(glue('{Sys.Date()} - Correlation matrix saved to {matrix_file}, Note that the correlation matrix might appear blank before it loads as this is a large file.'))
    
  }, error = function(e) {
    # Log and rethrow the error
    message(glue("ERROR: Failed to calculate correlation matrix: {e$message}"))
    stop(e)
  })
}



merge_CAGE_replicates <- function(myCAGEset, samples_tsv_file_path) {
  message(glue("{Sys.Date()} - Starting merging of replicates using {samples_tsv_file_path} ..."))
  samples <- tryCatch({
    read.delim(samples_tsv_file_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
  }, error = function(e) {
    stop(glue("Failed to read samples TSV file: {e$message}"))
  })
  
  if (!all(c("sample_name", "tissue", "genome") %in% colnames(samples))) {
    stop("ERROR: The samples TSV file must contain 'sample_name', 'tissue', and 'genome' columns.")
  }
  cage_sample_labels <- sampleLabels(myCAGEset)
  if (!all(samples$sample_name %in% cage_sample_labels)) {
    missing_samples <- setdiff(samples$sample_name, cage_sample_labels)
    stop(glue("ERROR: The following samples in the TSV file are missing in the CAGEset object: {paste(missing_samples, collapse = ', ')}"))
  }
  samples <- samples %>%
    mutate(mergedLabel = paste(tissue, genome, sep = "_"))
  mergeIndex <- match(samples$sample_name, cage_sample_labels)
  
  if (any(is.na(mergeIndex))) {
    stop("ERROR: Some CAGEset samples could not be matched to the metadata file. Check for inconsistencies.")
  }
  mergedSampleLabels <- unique(samples$mergedLabel)
  message(glue("Merging samples into groups: {paste(mergedSampleLabels, collapse = ', ')}"))
  
  myCAGEset <- tryCatch({
    mergeSamples(
      myCAGEset,
      mergeIndex = mergeIndex,
      mergedSampleLabels = samples$mergedLabel
    )
  }, error = function(e) {
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

normalise_and_cluster_samples <- function(myCAGEset, output_results_directory,  alpha, t, lower_range, upper_range, maxDist = 50, keepSingletonsAbove = 5) {
  # Normalises and the clusters the samples, returns the CAGEset object 
  # Saves the clustered 
  message(glue("{Sys.Date()} - Normalising samples with an alpha of {alpha}, T of {t}, lower_range of {lower_range} and upper_range of {upper_range} ..."))
  myCAGEset <- normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(lower_range, upper_range), alpha = alpha, T = t)
  message(glue("{Sys.Date()} - Normalising finished."))
  # Cluster
  message('Different clustering methods can be used, see word doc - \'CAGE_clustering_methods_lit_review\'')
  message(glue("{Sys.Date()} - Clustering samples with a maximum distance of {maxDist} and keeping singletons above {keepSingletonsAbove} ..."))
  myCAGEset <- distclu(myCAGEset, maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)
  message(glue("{Sys.Date()} - Clustering finished.") )
  # Save the clustered samples to a track
  # Save the CAGEset object with saveRDS
  saveRDS(myCAGEset, file = file.path(output_results_directory, "CAGE_object_clustered_samples.rds"))
  message(glue("{Sys.Date()} - Saved myCAGEset object to {output_results_directory}/CAGE_object_clustered_samples.rds with saveRDS."))
  message(glue("{Sys.Date()} - Exporting clustered samples to track ..."))
  smaple_names <- sampleLabels(myCAGEset)
  for (sample in smaple_names) {
    message(glue("{Sys.Date()} - Exporting {sample} ..."))
    sample_tagclusters <- tagClustersGR(myCAGEset, sample = sample)
    track_file_path <- file.path(output_results_directory, glue("{sample}_TC.bed"))
    trk <- exportToTrack(sample_tagclusters, what ="tagClusters",colorByExpressionProfile = FALSE,oneTrack = TRUE)
    export.bed(trk, con=track_file_path)
    message(glue("{Sys.Date()} - Exported  {sample} to track succsessfully."))
  }
  return(myCAGEset)
}

plot_location_tag_clusters <- function(myCAGEset, output_results_directory) {
  # Annotate Tag Clusters
  message(glue("{Sys.Date()} - Annotating tag clusters..."))
  genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
  myCAGEset <- annotateTagClusters(myCAGEset, genome)
  
  # Plot the annotation
  message(glue("{Sys.Date()} - Creating annotation plot..."))
  annot_plot <- plotAnnot(myCAGEset, what = "counts", main = "Distribution of CAGE Tags")
  plot_file_name <- file.path(output_results_directory, "plotAnnot_tag_clusters.pdf")
  ggsave(filename = plot_file_name, plot = annot_plot, width = 8, height = 6, dpi = 900)
  message(glue("{Sys.Date()} - Annotation plot saved to {plot_file_name}"))
}






main <- function() {
  args <- parse_arguments()
  myCAGEset <- read_in_cage_data(args$input_files, args$file_pattern) 
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
  # correlation_matrix - function works! just commenting it out while I work on other steps 
  correlation_matrix(myCAGEset, output_results_directory, tagCountThreshold = 5, applyThresholdBoth = FALSE, method="pearson") 

  myCAGEset <- merge_CAGE_replicates(myCAGEset, args$samples_tsv_file_path)
  # reverse_cumulatives - function works! just commenting it out while I work on other steps
  reverse_cumulatives(myCAGEset, output_results_directory, method = "powerLaw", fitInRange = c(5, 1000))
  saveRDS(myCAGEset, file ="myCAGEset_before_normalisation_and_clustering.rds")
  # Pause at this point to get the alpha and T, lower and upper values from the distribution plot
  # alpha = 1.06, T = 1e+07
  #alpha=1.06
  #t=1e+07
  #lower_range=1
  #upper_range=1000
  #myCAGEset <- normalise_and_cluster_samples(myCAGEset, output_results_directory, alpha, t, lower_range, upper_range, maxDist = 20, keepSingletonsAbove = 5)
  ## This works and you get an output file e.g. head /ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/R9_TC.bed
  #myCAGEset <- readRDS("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CAGE_object_clustered_samples.rds")
  #plot_location_tag_clusters_2(myCAGEset, output_results_directory)
 }


main()
traceback()
sessionInfo()

#saveRDS(myCAGEset, file ="testing.rds")
#read_in <- readRDS("testing.rds")