# The aim of this script is to determine if I can annotate TCs from a bed q file
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


#genome <- getBSgenome("BSgenome.Taestivum.ChineseSpring")
ce <- CAGEexp(genomeName = "BSgenome.Taestivum.ChineseSpring",
                     inputFiles = c('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/copy_no_merge_ALL_TC.bed'), 
                     inputFilesType = "bed",
                     sampleLabels = c('ALL'))

print(ce)
colData(ce)
ce <- getCTSS(ce)
CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
sampleLabels(ce)
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
ce <- annotateCTSS(ce, granges, upstream=500, downstream=500)
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
CTSScoordinatesGR(ce)
sample_tagcluster_gr <- CTSScoordinatesGR(ce)
# columns are "seqnames","pos","strand","genes","annotation" 

# Custom nearest-gene promoter function
nearest_gene_promoter <- function(seqnames, pos, strand, granges, distance = 500) {
  tc_gr <- GRanges(seqnames = Rle(seqnames),
                   ranges = IRanges(start = as.numeric(pos), width = 1),
                   strand = Rle(strand))
  nearest_hit <- nearest(tc_gr, granges)
  if (is.na(nearest_hit)) return("")
  nearest_gene_name <- as.character(mcols(granges)[nearest_hit, "gene_name"])
  nearest_gene_strand <- as.character(strand(granges)[nearest_hit])
  if (abs(pos - start(granges)[nearest_hit]) <= distance) return(nearest_gene_name)
  return("")
}

df <- as.data.frame(sample_tagcluster_gr)
head(df)
# Pre-filter for promoter TCs
promoter_rows <- df$annotation == "promoter"
promoter_df <- df[promoter_rows, ]
print(glue("Promoter TCs: {nrow(promoter_df)}"))

promoter_df$nearest_gene <- Map(function(seqnames, pos, strand) {
    nearest_gene_promoter(
      seqnames = as.character(seqnames),
      pos = as.numeric(pos),
      strand = as.character(strand),
      granges = granges_transcripts
    )
  }, promoter_df$seqnames, promoter_df$pos, promoter_df$strand)


df$nearest_gene <- ""
df$nearest_gene[promoter_rows] <- promoter_df$nearest_gene
original_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
df$gene <- ifelse(df$gene == "" & df$annotation == "promoter", df$nearest_gene, df$gene)
new_empty_genes <- sum(df$annotation == "promoter" & df$gene == "")
newly_annotated_genes <- original_empty_genes - new_empty_genes
print(glue("Added annotations to {newly_annotated_genes} promoter TCs."))
filename <- "test_my_annot_2.csv"
write.csv(df, filename, row.names = FALSE)
cat("Saved:", filename, "\n")


