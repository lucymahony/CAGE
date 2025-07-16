suppressPackageStartupMessages({
  library(ggplot2)
  library(topGO)
  library(dplyr)
  library(stringr)
})
#
## ---- CONFIG ----
input_dir <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data"
#gene2go_file <- file.path(input_dir, "geneID2GOrefseqv2.1.csv")
#background_file <- file.path(input_dir, "SIS_SLE_SRO_SSP_genes_shared_all.txt")
#files <- list.files(input_dir, pattern = "SIS_SLE_SRO_SSP_genes_unique_", full.names = TRUE)
#
## ---- LOAD GENE2GO ----
#geneID2GO <- readMappings(gene2go_file)
#geneID2GO <- lapply(geneID2GO, function(x) na.omit(x[x != ""]))
#geneUniverse_all <- names(geneID2GO)
#
## ---- GO ENRICHMENT FUNCTION ----
#run_topgo <- function(focus_genes, background_genes, geneID2GO, output_file, ontology = "BP") {
#  geneUniverse <- union(focus_genes, background_genes)
#  geneList <- factor(as.integer(geneUniverse %in% focus_genes))
#  table(geneList)
#  names(geneList) <- geneUniverse
#
#  if (length(unique(geneList)) < 2) {
#    message("Skipping ", output_file, ": geneList has only one class (not enough contrast)")
#    return(NULL)
#  }
#  GOdata <- new("topGOdata",
#                ontology = ontology,
#                allGenes = geneList,
#                nodeSize = 10,
#                annot = annFUN.gene2GO,
#                gene2GO = geneID2GO)
#
#  resultFisher <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
#  resultTable <- GenTable(GOdata,
#                          classicFisher = resultFisher,
#                          orderBy = "classicFisher",
#                          topNodes = min(200, length(score(resultFisher))))
#  resultTable$classicFisher <- as.numeric(gsub("< ", "", resultTable$classicFisher))
#  resultTable <- resultTable[resultTable$classicFisher <= 0.05, ]
#
#  write.csv(resultTable, file = output_file, row.names = FALSE)
#  message("GO terms saved to: ", output_file)
#}
#
## ---- RUN GO ENRICHMENT ----
#background_genes <- scan(background_file, what = "", sep = "\n")
#background_genes <- sub("\\.\\d+$", "", background_genes)  # Remove transcript version suffix
#background_genes <- background_genes[background_genes %in% geneUniverse_all]
#
## Create an empty list to collect results
#go_all_tissues <- list()
#
#for (file in files) {
#  cat("Processing:", file, "\n")
#  focus_genes <- scan(file, what = "", sep = "\n")
#  focus_genes <- sub("\\.\\d+$", "", focus_genes)
#  focus_genes <- focus_genes[focus_genes %in% geneUniverse_all]
#  
#  if (length(focus_genes) < 10) {
#    cat("Skipping ", basename(file), " (too few genes after filtering)\n")
#    next
#  }
#
#  tissue <- str_extract(basename(file), "(?<=unique_)[A-Za-z]+") %>% toupper()
#
#  # Collect results from all ontologies for this tissue
#  all_go <- lapply(c("BP", "MF", "CC"), function(ont) {
#    message("  Running GO enrichment for ontology: ", ont)
#    output_file_temp <- tempfile(fileext = ".csv")
#    run_topgo(focus_genes, background_genes, geneID2GO, output_file_temp, ontology = ont)
#    if (file.exists(output_file_temp)) {
#      df <- read.csv(output_file_temp)
#      df$ontology <- ont
#      return(df)
#    } else {
#      return(NULL)
#    }
#  }) %>% bind_rows()
#
#  # Tag with tissue name
#  if (nrow(all_go) > 0) {
#    all_go$tissue <- tissue
#    go_all_tissues[[tissue]] <- all_go
#  }
#}
#
## ---- After loop: Combine and write GO_ALL.csv ----
#combined_go <- bind_rows(go_all_tissues)
#if (nrow(combined_go) > 0) {
#  output_file <- file.path(input_dir, "GO_ALL.csv")
#  write.csv(combined_go, output_file, row.names = FALSE)
#  message("Saved full combined GO file: ", output_file)
#} else {
#  message("No GO results to write.")
#}
#
# ---- Bubble PLOT ----
# Define pattern to match your GO output files
go_file <- list.files(input_dir, pattern = "GO_ALL.csv", full.names = TRUE)
go_data <- read.csv(go_file)
print(head(go_data))
# Filter for ontology "BP"
go_data <- go_data %>%
  filter(ontology == "BP") 
# Filter and add -log2(p)
go_data <- go_data %>%
  filter(Significant >= 10, !is.na(classicFisher)) %>%
  mutate(logP = -log2(classicFisher))

# Sort by descening logP value -Take the smallest 20 terms e.g. top 20 rows of sorted_df
sorted_df <- go_data[order(-go_data$logP), ]
# Print the shape of the df
print('There are ')
print(dim(sorted_df))
print('identified terms')

top_terms <- head(sorted_df, 100)

# Bubble plot, faceted by ontology
p <- ggplot(top_terms, aes(x = tissue, y = Term)) +
  geom_point(aes(size = Significant, fill = logP),
             colour = "black", shape = 21, stroke = 0.25) +
  facet_wrap(~ontology, scales = "free_y") +
  scale_size_continuous(name = "Sig. genes") +
  scale_fill_gradient(low = "white", high = "#00A087FF", name = "-log2(p)") +
  xlab("Tissue") + ylab("GO term") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) + labs(title = "") + 
  scale_x_discrete(labels=c("SIS" = "Immature Spike", "SLE" = "Leaf", "SRO" = "Root", "SSP" = "Spike")) 

ggsave("GO_bubble_plot_BP.pdf", plot = p, width = 10, height = 12)
