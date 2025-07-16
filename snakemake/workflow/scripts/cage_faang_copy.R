# Code developed from: https://bitbucket.org/msalavat/cagewrap_public/src/master/CAGEfightR_aio.R


suppressPackageStartupMessages(library(CAGEfightR, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, quietly=T))
suppressPackageStartupMessages(library(BSgenome, quietly=T))
suppressPackageStartupMessages(library(BSgenome.Taestivum.ChineseSpring, quietly=T))
suppressPackageStartupMessages(library(tidyverse, quietly=T))
suppressPackageStartupMessages(library(magrittr, quietly=T))
suppressPackageStartupMessages(library(MultiAssayExperiment, quietly=T))
suppressPackageStartupMessages(library(SummarizedExperiment, quietly=T))
suppressPackageStartupMessages(library(Biostrings, quietly=T))
suppressPackageStartupMessages(library(GenomicFeatures, quietly=T))
suppressPackageStartupMessages(library(BiocParallel, quietly=T))
suppressPackageStartupMessages(library(InteractionSet, quietly=T))
suppressPackageStartupMessages(library(Gviz, quietly=T))
suppressPackageStartupMessages(library(tidyverse, quietly=T))
suppressPackageStartupMessages(library(FactoMineR, quietly=T))
suppressPackageStartupMessages(library(factoextra, quietly=T))
suppressPackageStartupMessages(library(reshape2, quietly=T))
suppressPackageStartupMessages(library(circlize, quietly=T))
suppressPackageStartupMessages(library(RColorBrewer, quietly=T))
suppressPackageStartupMessages(library(randomcoloR, quietly=T))
suppressPackageStartupMessages(library(philentropy, quietly=T))
suppressPackageStartupMessages(library(bioDist, quietly=T))
suppressPackageStartupMessages(library(patchwork, quietly=T))
suppressPackageStartupMessages(library(tidyverse, quietly=T))
suppressPackageStartupMessages(library(scales, quietly=T))
suppressPackageStartupMessages(library(tximport, quietly=T))
suppressPackageStartupMessages(library(Repitools, quietly=T))
suppressPackageStartupMessages(library(Hmisc, quietly=T))
suppressPackageStartupMessages(library(viridis, quietly=T))
suppressPackageStartupMessages(library(lme4, quietly=T))
suppressPackageStartupMessages(library(lmerTest, quietly=T))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(txdbmaker))

convert_gff3_to_txdb <- function(gff_file, output_file) {
    #' Convert GFF3 file to TxDb object
    #'
    #' This function converts a GFF3 file to a TxDb object, assigns transcript_id,
    #' removes exon phases, and saves the TxDb object to a specified file.
    #' For some reason some exons have a 0 in their phase rather than . 
    #' awk -F'\t' '$3 == "exon" && $8 == "."' $gff3 | wc -l
    #' 741975
    #' awk -F'\t' '$3 == "exon" && $8 == "0"' $gff3 | wc -l
    #' 250
    #' @param gff_file Path to the input GFF3 file.
    #' @param output_file Path to the output file where the TxDb object will be saved.
    #'
    #' @return None. The function saves the TxDb object to the specified file.
    gr <- import(gff_file, format="gff3")
    df <- as.data.frame(gr)
    df <- df %>%
      mutate(transcript_id = ifelse(type %in% c("mRNA", "transcript"), ID, NA))
    df <- df %>%
      group_by(Parent) %>%
      mutate(transcript_id = ifelse(is.na(transcript_id), first(na.omit(transcript_id)), transcript_id)) %>%
      ungroup()
    mcols(gr)$transcript_id <- df$transcript_id
    mcols(gr)$phase <- ifelse(gr$type == "exon", NA, mcols(gr)$phase)
    txdb <- makeTxDbFromGRanges(gr)
    saveDb(txdb, file=output_file)
    }

gff_file <- '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3'
txdb_file <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.txdb"
#convert_gff3_to_txdb(gff_file, txdb_file)
# Read in the TxDb object
txdb <- loadDb(txdb_file)


setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/")
bw_plus <- BigWigFileList(c("CL1.CS.se.star_plus.unique.bw", "CL2.CS.se.star_plus.unique.bw", "CL3.CS.se.star_plus.unique.bw", "CL4.CS.se.star_plus.unique.bw", "CL5.CS.se.star_plus.unique.bw"))
bw_minus <- BigWigFileList(c("CL1.CS.se.star_minus.unique.bw", "CL2.CS.se.star_minus.unique.bw", "CL3.CS.se.star_minus.unique.bw", "CL4.CS.se.star_minus.unique.bw", "CL5.CS.se.star_minus.unique.bw"))
names(bw_minus) <- c("CL1", "CL2", "CL3", "CL4", "CL5")
names(bw_plus) <- c("CL1", "CL2", "CL3", "CL4", "CL5")

# Extract and print sequence names and lengths for each BigWigFile

print(seqnames(seqinfo(bw_plus)))
print(seqlengths(seqinfo(bw_plus)))

print(seqnames(seqinfo(bw_minus)))
print(seqlengths(seqinfo(bw_minus)))


##As the conversion of the bigWigs are done using 1 indexed coordinates the legth of the BSgenome is 1 bp short of the bw files. Fixed it manually.
chinese_spring_seqinfo <- seqinfo(BSgenome.Taestivum.ChineseSpring)
BSgenome.Taestivum.ChineseSpring@seqinfo@seqlengths <- as.integer(seqlengths(chinese_spring_seqinfo) + 1)
print(seqlengths(BSgenome.Taestivum.ChineseSpring))

print('check here: ')
print(seqnames(seqinfo(bw_plus)))
print(seqnames(seqinfo(BSgenome.Taestivum.ChineseSpring)))
Gm <- seqinfo(bw_plus[[1]])

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome = Gm) 


print(CTSSs)

CTSSs<-CTSSs %>% calcTPM() %>% calcPooled()
TCs<-quickTSSs(CTSSs)

TSSs <- TCs %>%
  calcTPM() %>%
  subsetBySupport(inputAssay="counts",
                  unexpressed=1, #minimum TPM requirement
                  minSamples=2) %>%  #The 2/3rd representation criteria. 
  calcTPM()
print(TSSs)


TSSs_full <- TCs %>%
  calcTPM() %>%
  subsetBySupport(inputAssay="counts",
                  unexpressed=1, #minimum TPM requirement
                  minSamples=1) %>% #The tissue specific set. 
  calcTPM()

print(TSSs_full)


# Annotation 

TSSs <- assignTxType(TSSs,txModels = txdb,swap="thick")
TSSs<-assignTxID(object = TSSs,txModels = txdb,swap="thick")
TSSs <- assignGeneID(TSSs, txdb)
TSSs <- assignMissingID(TSSs, outputColumn = "geneID", prefix = "NOVELG")
TSSs <- assignMissingID(TSSs, outputColumn = "txID", prefix = "NOVELT")


TSSs_df<-TSSs %>% rowRanges() %>% data.frame()

##saving point
con <- pipe("xz -T8 -6 -e > crf_TCs_object.xz", "wb")
save(TCs, file = con); close(con)
rm(con)

con <- pipe("xz -T8 -6 -e > crf_TSSs_object.xz", "wb")
save(TSSs, file = con); close(con)
rm(con)

con <- pipe("xz -T8 -6 -e > crf_CTSS_TPM_object.xz", "wb")
save(CTSSs, file = con); close(con)
rm(con)



#CTSSs<-CTSSs %>% calcTPM() %>% calcPooled()
#TCs<-quickTSSs(CTSSs)
#
#TSSs <- TCs %>%
#  calcTPM() %>%
#  subsetBySupport(inputAssay="counts",
#                  unexpressed=10, #minimum TPM requirement
#                  minSamples=2) %>%  #The 2/3rd representation criteria. 
#  calcTPM()
#
#
#TSSs_full <- TCs %>%
#  calcTPM() %>%
#  subsetBySupport(inputAssay="counts",
#                  unexpressed=10, #minimum TPM requirement
#                  minSamples=1) %>% #The tissue specific set. 
#  calcTPM()
#

#read_bigwig_files <- function(path, pattern) {
#  files <- dir(path = path, pattern = pattern, full.names = TRUE, include.dirs = TRUE)
#  return(BigWigFileList(files))
#}
#
#save_object <- function(obj, filename) {
#  con <- pipe(paste0("xz -T8 -6 -e > ", filename, ".xz"), "wb")
#  save(obj, file = con)
#  close(con)
#}
#
#perform_pca <- function(data, n_comp = 5) {
#  PCA(data, scale.unit = TRUE, ncp = n_comp, graph = FALSE)
#}
#
## Load tally file
#tally <- read.csv(config$tally_file, header = TRUE, stringsAsFactors = FALSE) %>%
#  clean_tissue_names()
#
## Process BAM files
#bam_files <- process_bam_files(config$input_dir)
#bw_plus <- read_bigwig_files(config$bw_dir, ".plus.bw")
#bw_minus <- read_bigwig_files(config$bw_dir, ".minus.bw")
#
#CTSSs <- quantifyCTSSs(
#  plusStrand = bw_plus,
#  minusStrand = bw_minus,
#  genome = SeqinfoForBSGenome(BSgenome.Oaries.NCBI.Ramb1.0)
#) %>%
#  calcTPM() %>%
#  calcPooled()
#
#TSSs_PCA <- perform_pca(TSSs_pca_input)
#BCs_PCA <- perform_pca(BCs_pca_input)
#
#save_object(TCs, "crf_TCs_object")
#save_object(TSSs, "crf_TSSs_object")
#save_object(CTSSs, "crf_CTSS_TPM_object")
#
#cage_write <- function(df, tissue, filename_prefix) {
#  tmp <- data.frame(
#    cluster = row.names(df),
#    tissue = df[, tissue]
#  ) %>%
#    filter(tissue != 0)
#
#  tmp1 <- str_split_fixed(tmp$cluster, ":", 2)
#  tmp2 <- str_split_fixed(tmp1[, 2], ";", 2)
#  tmp3 <- str_split_fixed(tmp2[, 1], "-", 2)
#
#  ready <- data.frame(
#    chr = tmp1[, 1],
#    start = tmp3[, 1],
#    end = tmp3[, 2],
#    strand = tmp2[, 2],
#    TPM = tmp$tissue
#  )
#
#  write.table(
#    ready,
#    paste0(config$output_dir, "/", filename_prefix, "_", tissue, "_tpm.tsv"),
#    quote = FALSE,
#    sep = "\t",
#    row.names = FALSE
#  )
#}
#
## Apply to all tissues
#walk(colnames(TSSs_pca_input), ~cage_write(TSSs_pca_input, .x, "TSSs"))
#walk(colnames(BCs_pca_input), ~cage_write(BCs_pca_input, .x, "BCs"))
#
#calculate_mi <- function(data) {
#  data_invert <- t(data)
#  colnames(data_invert) <- data$txID
#  return(MIdist(data_invert))
#}
#
#RNA_MI <- calculate_mi(RNA)
#TSSs_MI <- calculate_mi(TSSs_pca_input)
#BCs_MI <- calculate_mi(BCs_pca_input)
#
#plot_dendrogram <- function(mi_matrix, title) {
#  factoextra::fviz_dend(
#    hclust(mi_matrix, "complete"),
#    k = 12,
#    k_colors = "lancet",
#    type = "phylogenic",
#    repel = TRUE,
#    color_labels_by_k = TRUE,
#    cex = 0.7,
#    phylo_layout = "layout.auto"
#  ) + ggtitle(title)
#}
#
#cowplot::plot_grid(
#  plot_dendrogram(RNA_MI, "RNA-Seq TPM Clusters"),
#  plot_dendrogram(TSSs_MI, "Unidirectional TSSs Clusters"),
#  labels = c("RNA-Seq TPM Clusters", "Unidirectional TSSs Clusters")
#)
#