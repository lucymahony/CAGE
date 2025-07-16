sink(file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/CS.star.CAGE_fight_r_analysis.txt")

library(GenomicAlignments)
library(rtracklayer)
library(CAGEfightR)
library(glue)
library(GenomicFeatures)
library(dplyr)


parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4) {
    stop("ERROR: Please provide four commandline arguments: <bw_plus_files> <bw_minus_files> <input_files_directory> <samples_tsv_file_path>")
  }
  bw_plus_files_pattern <- args[1]
  bw_minus_files_pattern <- args[2] 
  input_files_directory <- args[3]
  samples_tsv_file_path <- args[4]
  message("WARNING: The successful running of this script requires the config files to be correctly formatted.")
  message("Check the samples TSV file contains the samples you wish to process ")
  list(
    input_files_directory = input_files_directory,
    samples_tsv_file_path = samples_tsv_file_path
  )
}

validate_and_clean_gff3 <- function(gff_file, output_file = NULL) {
  gff <- import(gff_file, format = "gff3")
  features_to_check <- c("gene", "mRNA", "exon", "CDS")
  missing_transcript_id <- list()
  for (feature in features_to_check) {
    subset <- gff[gff$type == feature]
    missing <- subset[is.na(subset$transcript_id)]
    if (length(missing) > 0) {
      message("Features of type '", feature, "' missing 'transcript_id':")
      print("Number of missing features:")
      print(length(missing))
      print("Details of missing features:")
      print(missing)
      if ("ID" %in% names(mcols(missing))) {
        message("IDs of missing features:")
        print(missing$ID)
      } else {
        message("No 'ID' attribute found in missing features.")
      }
    } else {
      message("No features of type '", feature, "' are missing 'transcript_id'.")
    }
  }
  if (all(sapply(features_to_check, function(feature) {
    subset <- gff[gff$type == feature]
    sum(is.na(subset$transcript_id)) == 0
  }))) {
    message("All features have 'transcript_id' attributes.")
  }
  transcript_ids <- gff$transcript_id[!is.na(gff$transcript_id)]
  duplicate_ids <- transcript_ids[duplicated(transcript_ids)]
  if (length(duplicate_ids) > 0) {
    warning("Duplicate 'transcript_id' values found:")
    print(unique(duplicate_ids))
  } else {
    message("All 'transcript_id' values are unique.")
  }
  exon_indices <- which(gff$type == "exon")
  if (any(!is.na(gff$phase[exon_indices]))) {
    warning("'phase' attribute found in exon features. Removing...")
    gff$phase[exon_indices] <- NA
  }
  if (!is.null(output_file)) {
    export(gff, output_file, format = "gff3")
  }
  return(gff)
}

create_TxDb <- function(gff_file) {
  txdb <- makeTxDbFromGFF(gff_file)
  return(txdb)
}

quantifying_and_annotating_tss_and_enhancers <- function(txdb) {
    setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/copy_of_CS_se_star_sorted_bw")
    bw_plus <- BigWigFileList(c("CL1.CS.se.star_plus.unique.bw","CL2.CS.se.star_plus.unique.bw", "CL3.CS.se.star_plus.unique.bw", "CL4.CS.se.star_plus.unique.bw", "CL5.CS.se.star_plus.unique.bw",
                         "CR1.CS.se.star_plus.unique.bw", "CR2.CS.se.star_plus.unique.bw", "CR3.CS.se.star_plus.unique.bw", "CR4.CS.se.star_plus.unique.bw", "CR5.CS.se.star_plus.unique.bw", 
                        "SIS1.CS.se.star_plus.unique.bw", "SIS2.CS.se.star_plus.unique.bw", "SIS3.CS.se.star_plus.unique.bw",
                         "SLE1.CS.se.star_plus.unique.bw", "SLE2.CS.se.star_plus.unique.bw", "SLE3.CS.se.star_plus.unique.bw",
                          "SRO1.CS.se.star_plus.unique.bw", "SRO2.CS.se.star_plus.unique.bw", "SRO3.CS.se.star_plus.unique.bw",
                           "SSP1.CS.se.star_plus.unique.bw", "SSP2.CS.se.star_plus.unique.bw", "SSP3.CS.se.star_plus.unique.bw"))
    bw_minus <- BigWigFileList(c("CL1.CS.se.star_minus.unique.bw","CL2.CS.se.star_minus.unique.bw", "CL3.CS.se.star_minus.unique.bw", "CL4.CS.se.star_minus.unique.bw", "CL5.CS.se.star_minus.unique.bw",
                                "CR1.CS.se.star_minus.unique.bw", "CR2.CS.se.star_minus.unique.bw", "CR3.CS.se.star_minus.unique.bw", "CR4.CS.se.star_minus.unique.bw", "CR5.CS.se.star_minus.unique.bw",
                                "SIS1.CS.se.star_minus.unique.bw", "SIS2.CS.se.star_minus.unique.bw", "SIS3.CS.se.star_minus.unique.bw",
                                "SLE1.CS.se.star_minus.unique.bw", "SLE2.CS.se.star_minus.unique.bw", "SLE3.CS.se.star_minus.unique.bw",
                                "SRO1.CS.se.star_minus.unique.bw", "SRO2.CS.se.star_minus.unique.bw", "SRO3.CS.se.star_minus.unique.bw",
                                "SSP1.CS.se.star_minus.unique.bw", "SSP2.CS.se.star_minus.unique.bw", "SSP3.CS.se.star_minus.unique.bw"))
    names(bw_plus) <- c("CL1", "CL2", "CL3", "CL4", "CL5", "CR1", "CR2", "CR3", "CR4", "CR5", "SIS1", "SIS2", "SIS3", "SLE1", "SLE2", "SLE3", "SRO1", "SRO2", "SRO3", "SSP1", "SSP2", "SSP3")
    names(bw_minus) <- c("CL1", "CL2", "CL3", "CL4", "CL5", "CR1", "CR2", "CR3", "CR4", "CR5", "SIS1", "SIS2", "SIS3", "SLE1", "SLE2", "SLE3", "SRO1", "SRO2", "SRO3", "SSP1", "SSP2", "SSP3")
    Gm <- seqinfo(bw_plus[[1]])
    CTSSs <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, genome=Gm)
    print(CTSSs)
    rowRanges(CTSSs)
    CTSSs <- CTSSs %>%
    calcTPM() %>%
    calcPooled()
    TCs <- quickTSSs(CTSSs)
    print(TCs)
    TSSs <- TCs
    print('No subsetting by support for now')
    #TSSs <- TCs %>%
    #  calcTPM() %>%
    #  subsetBySupport(inputAssay = "TPM", 
    #                  unexpressed = 1, 
    #                  minSamples = 2)
    print(TSSs)
    BCs <- quickEnhancers(CTSSs)
    
    TSSs <- assignTxID(TSSs, txModels = txdb, swap = "thick")
    TSSs <- assignTxType(TSSs, txModels = txdb, swap = "thick")
    BCs <- assignTxType(BCs, txModels = txdb, swap = "thick")
    Enhancers <- subset(BCs, txType %in% c("intergenic", "intron"))

    # merging
    TSSs$totalTags <- NULL
    Enhancers$totalTags <- NULL
    rowData(TSSs)$balance <- NA
    rowData(TSSs)$bidirectionality <- NA
    rowData(Enhancers)$txID <- NA

    rowData(TSSs)$clusterType <- "TSS"
    rowData(Enhancers)$clusterType <- "Enhancer"
    RSE <- combineClusters(object1 = TSSs, 
                       object2 = Enhancers, 
                       removeIfOverlapping = "object1")
    RSE <- calcTPM(RSE)
    axis_track <- GenomeAxisTrack()
    tx_track <- GeneRegionTrack(txdb, 
                                name = "Gene Models", 
                                col = NA,
                                fill = "bisque4", 
                                shape = "arrow", 
                                showId = TRUE)
    plot_region <- RSE %>% 
    rowRanges() %>% 
    subset(clusterType == "TSS") %>% 
    .[1] %>%
    add(100) %>%
    unstrand()
    ctss_track <- CTSSs %>%
        rowRanges() %>%
        subsetByOverlaps(plot_region) %>%
        trackCTSS(name = "CTSSs")
    cluster_track <- RSE %>%
        subsetByOverlaps(plot_region) %>%
        trackClusters(name = "Clusters", 
                      col = NA, 
                      showId = TRUE)
    plotTracks(list(axis_track, 
                    ctss_track,
                    cluster_track,
                    tx_track),
               from = start(plot_region), 
               to = end(plot_region), 
               chromosome = as.character(seqnames(plot_region)))

    # save the plot
    ggsave("TSSs.pdf", width = 8, height = 6)

    plot_region <- RSE %>% 
        rowRanges() %>% 
        subset(clusterType == "Enhancer") %>% 
        .[1] %>%
        add(100) %>%
        unstrand()
    ctss_track <- CTSSs %>%
        rowRanges() %>%
        subsetByOverlaps(plot_region) %>%
        trackCTSS(name = "CTSSs")
    cluster_track <- RSE %>%
        rowRanges %>%
        subsetByOverlaps(plot_region) %>%
        trackClusters(name = "Clusters", 
                      col = NA, 
                      showId = TRUE)
    plotTracks(list(axis_track, 
                    ctss_track,
                    cluster_track,
                    tx_track),
               from = start(plot_region), 
               to = end(plot_region), 
               chromosome = as.character(seqnames(plot_region)))
    # save the plot
    ggsave("Enhancers.pdf", width = 8, height = 6)

    cluster_info <- RSE %>%
      rowData() %>%
      as.data.frame()

    ggplot(cluster_info, aes(x = txType, fill = clusterType)) +
      geom_bar(alpha = 0.75, position = "dodge", color = "black") +
      scale_fill_colorblind("Cluster type") +
      labs(x = "Cluster annotation", y = "Frequency") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    # save the plot
    ggsave("Cluster_info.pdf", width = 8, height = 6)

    ggplot(cluster_info, aes(x = txType, 
                           y = log2(score / ncol(RSE)), 
                           fill = clusterType)) +
      geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.50, 0.75)) +
      scale_fill_colorblind("Cluster type") +
      labs(x = "Cluster annotation", y = "log2(TPM)") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave("Cluster_info_violin.pdf", width = 8, height = 6)

}

analysing_tss_shapes_and_sequences <- function(){

}

finding_interacting_enhances_and_tss <- function(){

}
finding_streches_for_enhancers <- function(){

}



main <- function() {
    #args <- parse_arguments()
    #print(args)
    #input_files_directory<- args$input_files_directory
    #samples_tsv_file_path <- args$samples_tsv_file_path
    #txdb <- create_TxDb()
    #print(txdb)
    #read_in_cage_data(txdb)
    gff_file <- '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3'
    #cleaned_gff <- 'iwgsc_refseqv2.1_annotation_200916_HC_cleaned.gff3'
    #validate_and_clean_gff3(gff_file, cleaned_gff)
    txdb <- create_TxDb(gff_file)
    quantifying_and_annotating_tss_and_enhancers(txdb)
    

}

main()
sink(file = NULL)