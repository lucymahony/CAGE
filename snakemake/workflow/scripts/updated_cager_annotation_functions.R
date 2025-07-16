# Custom Promoter Annotation Function 

# In the package CAGEr, there is an issue with annotating TC's in the promoter region, because of the function 
# ranges2genes doesn't consider up and downstream parameters, where as ranges2annot does so you end up with TCs being assigned 'promoter' with no gene name. 
# Run this script post annotateTagClusters to fix this issue.

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
})

setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data")
output_results_directory <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps"
#Distclu_myCAGEset <- readRDS("no_merge_reps/distclu_myCAGEset.rds")
#
## GFF Granges required 
#convert_gff_to_granges <- function(gff_file, transcript_type = "protein_coding") {
#  print('Generating own GRanges object from GFF file .... ')
#  data_frame <- read.delim(gff_file, header = FALSE, skip = 1, stringsAsFactors = FALSE, colClasses = c(rep("character", 8), "character"))
#  colnames(data_frame) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
#  data_frame$feature[data_frame$feature == "mRNA"] <- "transcript"
#  if (!is.character(data_frame$attributes)) {
#    data_frame$attributes <- as.character(data_frame$attributes)
#  }
#  data_frame <- data_frame %>%
#    mutate(
#      gene_name = case_when(
#        grepl("ID=([^;]+)", attributes) ~ sub(".*ID=([^;]+).*", "\\1", attributes),
#        grepl("gene_id=([^;]+)", attributes) ~ sub(".*gene_id=([^;]+).*", "\\1", attributes),
#        grepl("Parent=([^;]+)", attributes) ~ sub(".*Parent=([^;]+).*", "\\1", attributes),
#        TRUE ~ NA_character_  # Assign NA if no match is found
#      )
#    )
#  gr <- GRanges(
#    seqnames = Rle(data_frame$seqname),
#    ranges = IRanges(start = as.numeric(data_frame$start), end = as.numeric(data_frame$end)),
#    strand = Rle(data_frame$strand),
#    type = Rle(data_frame$feature),
#    transcript_type = Rle(rep(transcript_type, nrow(data_frame))),
#    gene_name = Rle(data_frame$gene_name)
#  )
#  return(gr)
#}
#granges <- convert_gff_to_granges("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3")
#
#distclu_myCAGEset <- annotateCTSS(distclu_myCAGEset, granges, upstream=500, downstream=500)
#distclu_myCAGEset <- annotateTagClusters(distclu_myCAGEset, granges, upstream = 500, downstream = 500)
#print(tagClustersGR(distclu_myCAGEset, 1)) # There should be 5 metadata columns score, dominant_ctss, nr_ctss, annotation, genes - there actually only seems to be  3 metadata columns
## Turn tagClustersGR(distclu_myCAGEset, 1) into a data frame
#tagClustersGR_df <- as.data.frame(tagClustersGR(distclu_myCAGEset, 1))
#print(glue("Number of rows with 'promoter' annotation: {nrow(tagClustersGR_df[tagClustersGR_df$annotation == 'promoter', ])}"))
#print(glue("Number of rows with 'promoter' annotation and genes = '': {nrow(tagClustersGR_df[tagClustersGR_df$annotation == 'promoter' & tagClustersGR_df$genes == '', ])}"))
#
#
## Define function to find nearest gene
#nearest_gene_promoter <- function(annotation, genes, seqnames, start, end, strand, dominant_ctss_pos, granges) {
#  if (annotation != "promoter") {
#    return("")
#  }
#  tc_gr <- GRanges(seqnames = Rle(seqnames),
#                   ranges = IRanges(start = as.numeric(start), end = as.numeric(end)),
#                   strand = Rle(strand))
#  granges <- granges[granges$type == "transcript"]
#  nearest_hit <- nearest(tc_gr, granges)
#  if (is.na(nearest_hit)) {
#    return("")
#  }
#  nearest_gene_name <- as.character(mcols(granges)[nearest_hit, "gene_name"])
#  nearest_gene_strand <- as.character(strand(granges)[nearest_hit])
#  if (strand != nearest_gene_strand) {
#    return("")
#  } else if (abs(dominant_ctss_pos - start(granges)[nearest_hit]) <= 500) {
#    return(nearest_gene_name)
#  } else {
#    return("")
#  }
#}
#
## Assin nearest gene 
#tagClustersGR_df$nearest_gene <- ""
#tagClustersGR_df$nearest_gene <- Map(function(annotation, genes, seqnames, start, end, strand, dominant_ctss_pos) {
#  nearest_gene_promoter(
#    annotation = as.character(annotation),
#    genes = as.character(genes),
#    seqnames = as.character(seqnames),
#    start = as.numeric(start),
#    end = as.numeric(end),
#    strand = as.character(strand),
#    dominant_ctss_pos = as.numeric(dominant_ctss_pos),
#    granges = granges
#  )
#}, tagClustersGR_df$annotation, tagClustersGR_df$genes, tagClustersGR_df$seqnames, 
#   tagClustersGR_df$start, tagClustersGR_df$end, tagClustersGR_df$strand, 
#   tagClustersGR_df$dominant_ctss.pos)
#
## Turn the unassigned data frame into a GRanges object
#tagClustersGR_df$seqnames <- as.character(tagClustersGR_df$seqnames)
#tagClustersGR_df$start <- as.numeric(tagClustersGR_df$start)
#tagClustersGR_df$end <- as.numeric(tagClustersGR_df$end)
#tagClustersGR_df$strand <- as.character(tagClustersGR_df$strand)
#tagClustersGR_df$score <- as.numeric(tagClustersGR_df$score)
#tagClustersGR_df$dominant_ctss_pos <- as.numeric(tagClustersGR_df$dominant_ctss.pos)
#tagClustersGR_df$dominant_ctss_score <- as.numeric(tagClustersGR_df$dominant_ctss.score)
#tagClustersGR_df$nr_ctss <- as.numeric(tagClustersGR_df$nr_ctss)
#tagClustersGR_df$annotation <- as.character(tagClustersGR_df$annotation)
#tagClustersGR_df$genes <- as.character(tagClustersGR_df$genes)
#tagClustersGR_df$nearest_gene <- as.character(tagClustersGR_df$nearest_gene)
#tagClustersGR_df$nearest_gene[is.na(tagClustersGR_df$nearest_gene)] <- "" # Turn Nas into ""
#
#tagClustersGR_gr <- GRanges(
#  seqnames = Rle(tagClustersGR_df$seqnames),
#  ranges = IRanges(start = tagClustersGR_df$start, end = tagClustersGR_df$end),
#  strand = Rle(tagClustersGR_df$strand),
#  score = Rle(tagClustersGR_df$score),
#  dominant_ctss_pos = Rle(tagClustersGR_df$dominant_ctss_pos),
#  dominant_ctss_score = Rle(tagClustersGR_df$dominant_ctss_score),
#  nr_ctss = Rle(tagClustersGR_df$nr_ctss),
#  annotation = Rle(tagClustersGR_df$annotation),
#  genes = Rle(tagClustersGR_df$genes),
#  nearest_gene = Rle(tagClustersGR_df$nearest_gene)
#  )
#
## Save the GRanges object
#print(head(tagClustersGR_gr))
#print('saving ....')
#saveRDS(tagClustersGR_gr, file = "no_merge_reps/unassigned_gr.rds")

# Check that the nearest gene has been assigned correctly
tagClustersGR_gr <- readRDS("no_merge_reps/unassigned_gr.rds")
df <- as.data.frame(tagClustersGR_gr)
original_empty_genes <- nrow(df[df$annotation == 'promoter' & df$gene == '', ])
newly_annotated_genes <- nrow(df[df$annotation == 'promoter' & df$gene == '' & df$nearest_gene !='', ])

# If the gene is empty then assign it to the nearest gene, when the annotation is promoter
print('Assigning nearest gene to promoter TCs with empty gene metadata column')
df$gene <- ifelse(df$gene == "" & df$annotation=='promoter', df$nearest_gene, df$gene)
# To check that the merging worked correctly, the number of empty promoter genes should be 
# The original empty - the number where the nearest gene is not empty but gene is.
new_empty_genes <- nrow(df[df$annotation == 'promoter' & df$gene == '', ])

print(glue("Original empty genes: {original_empty_genes}, added the annotation of {newly_annotated_genes} genes, which means there are now only {new_empty_genes} empty genes in the promoter annotation."))
# Check original_empty_genes - newly_annotated_genes == new_empty_genes

if (original_empty_genes - newly_annotated_genes == new_empty_genes) {
  print('Wooo!')
} else {
  print('Something went wrong!')
}
