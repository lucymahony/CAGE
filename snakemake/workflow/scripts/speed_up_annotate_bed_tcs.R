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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript script.R <BASE_DIR> <BED_FILENAME>")
}
base_dir <- args[1]
bed_file_name <- args[2]
# Paths
bed_file <- file.path(base_dir, bed_file_name)
sample_label <- sub("_.*", "", bed_file_name) # Sample label derived from BED file (e.g., "CL_TC.bed" → "CL")
output_file_path <- file.path(base_dir, paste0(sample_label, "_annotated.csv"))

ce <- CAGEexp(
  genomeName = "BSgenome.Taestivum.ChineseSpring",
  inputFiles = c(bed_file), 
  inputFilesType = "bed",
  sampleLabels = c(sample_label)
)

# Get CTSSs
ce <- getCTSS(ce)

# Import and filter gene annotations (faster and cleaner)
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

# Annotate CTSSs
ce <- annotateCTSS(ce, granges, upstream = 500, downstream = 500)

# Extract CTSS info
tc_gr <- CTSScoordinatesGR(ce)
df <- as.data.frame(tc_gr)

# Identify promoter rows
promoter_rows <- df$annotation == "promoter"
promoter_df <- df[promoter_rows, ]
cat(glue("Promoter TCs: {nrow(promoter_df)}\n"))

# Vectorized nearest transcript search for promoter-tagged CTSS
promoter_tc_gr <- GRanges(
  seqnames = promoter_df$seqnames,
  ranges = IRanges(start = as.numeric(promoter_df$pos), width = 1),
  strand = promoter_df$strand
)

nearest_hits <- nearest(promoter_tc_gr, granges_transcripts)
nearest_gene_names <- rep("", length(promoter_tc_gr))
valid_hits <- !is.na(nearest_hits)
nearest_gene_names[valid_hits] <- mcols(granges_transcripts)$gene_name[nearest_hits[valid_hits]]

# Merge back into dataframe
df$nearest_gene <- ""
df$nearest_gene[promoter_rows] <- nearest_gene_names

# Fill in missing gene annotations for promoters
original_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
df$gene <- ifelse(df$gene == "" & df$annotation == "promoter", df$nearest_gene, df$gene)
new_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
newly_annotated_genes <- original_empty_genes - new_empty_genes

cat(glue("Added annotations to {newly_annotated_genes} promoter TCs.\n"))

# Save results
write.csv(df, output_file_path, row.names = FALSE)
cat("Saved:", output_file_path, "\n")
