
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

annotating_lc_genes <- function(list_input_files, gff_file_lc, upstream, downstream, output_results_directory) {
  # For prior generated TagClusters, they are read in using existing function read_in_cage_data pretending the tag clusters are 'ctss'
  # And saves the annotated LC genes to a csv file for each tissue/genome TagCluster file
  granges <- convert_gff_to_granges(gff_file_lc, transcript_type = "protein_coding")
  for (file in list_input_files) {
    message(glue("{Sys.Date()} - Annotating LC genes for file {file}"))
    file_pattern <-"_at_least_2_overlaps_max_score.bed"
    myCAGEset <- read_in_cage_data(file, file_pattern)
    # Name the csv file after the input files
    csv_file_name <- sub(file_pattern, "", basename(file))
    csv_file_name <- paste0(csv_file_name, "_annotated_lc_genes_ctss.csv")
    annotated_ctss(myCAGEset, granges, output_results_directory, csv_file_name, upstream, downstream)
  }
}

count_lc_genes_with_tag_clusters <- function(list_input_files, output_results_directory, upstream, downstream) {
  # Uses the annotated LC genes csv files to count the number of LC genes that are supported by TC's.
  gene_tissue_map <- list()
  for (file in list_input_files) {
    tissue_name <- sub("_at_least_2_overlaps_max_score.bed", "", basename(file))
    message(glue("{Sys.Date()} - Counting LC genes with tag clusters for file {file}"))
    csv_file_name <- paste0("up_", upstream, "_down_", downstream, "_", tissue_name, "_annotated_lc_genes_ctss.csv")
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

  ## List the input files as all the files in the output directory that end in   at_least_2_overlaps_max_score.bed
  list_input_files <- list.files("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data", pattern = "_at_least_2_overlaps_max_score.bed", full.names = TRUE)
  #gff_file_lc <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_LC.gff3"
  output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/"
  ## Generate the annotated CTSS data
  #annotating_lc_genes(list_input_files, gff_file_lc, 500, 500, output_results_directory)
  ## Read in the annotated CTSS data, and subset csv annotation = promoter, list the ‘genes’ column, get unique gene list and count the number of unique genes

  gene_tissue_df <- count_lc_genes_with_tag_clusters(list_input_files, output_results_directory, 500, 500)


}

main()

