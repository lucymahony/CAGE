# Load libraries
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
  library(data.table)
})

setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data")
myCAGEset <- readRDS("distclu_myCAGEset.rds")
CTSS_gr <- CTSScoordinatesGR(myCAGEset)


# Extract promoter regions from CAGEexp object
# Define a list of sample names
samples <- c("SRO1", "SRO2", "SRO3", "SLE1", "SLE2", "SLE3", "SSP1", "SSP2", "SSP3", "SIS1", "SIS2", "SIS3", "CR1", "CR2", "CR3", "CR4", "CR5", "CL1", "CL2", "CL3", "CL4", "CL5")
file_path <- "promoter_comparison_from_distclu_myCAGEset.csv"

# Extract promoter regions for each sample
promoter_list <- lapply(samples, function(sample) tagClustersGR(myCAGEset, sample = sample))
names(promoter_list) <- samples

# Compute pairwise intersections and Jaccard similarity
pairwise_shared <- combn(samples, 2, function(pair) {
  sample1 <- pair[1]
  sample2 <- pair[2]
  
  intersect_size <- length(intersect(promoter_list[[sample1]], promoter_list[[sample2]]))
  union_size <- length(promoter_list[[sample1]]) + length(promoter_list[[sample2]]) - intersect_size
  jaccard <- ifelse(union_size > 0, intersect_size / union_size, NA)
  
  list(shared = intersect_size, jaccard = jaccard)
}, simplify = FALSE)
names(pairwise_shared) <- apply(combn(samples, 2), 2, paste, collapse = " and ")

# Compute shared across all samples
shared_all <- length(Reduce(intersect, promoter_list))

# Compute unique promoters for each sample
unique_counts <- sapply(samples, function(sample) {
  others <- setdiff(samples, sample)
  length(setdiff(promoter_list[[sample]], Reduce(union, promoter_list[others])))
})
names(unique_counts) <- paste("Unique to", samples)

# Create a data frame for CSV output
upset_data <- data.frame(
  Sample1 = character(),
  Sample2 = character(),
  Shared = integer(),
  Jaccard = numeric()
)

for (name in names(pairwise_shared)) {
  pair <- unlist(strsplit(name, " and "))
  upset_data <- rbind(upset_data, data.frame(
    Sample1 = pair[1],
    Sample2 = pair[2],
    Shared = pairwise_shared[[name]]$shared,
    Jaccard = pairwise_shared[[name]]$jaccard
  ))
}

# Add unique counts for individual samples
unique_data <- data.frame(
  Sample1 = names(unique_counts),
  Sample2 = "Unique",
  Shared = unique_counts,
  Jaccard = NA
)

# Combine and save to CSV
final_data <- rbind(upset_data, unique_data)
write.csv(final_data, file_path, row.names = FALSE)

# Print results
for (name in names(pairwise_shared)) {
  cat(name, "Shared:", pairwise_shared[[name]]$shared, "Jaccard Similarity:", pairwise_shared[[name]]$jaccard, "\n")
}
cat("Shared across all samples:", shared_all, "\n")
for (name in names(unique_counts)) {
  cat(name, ":", unique_counts[[name]], "\n")
}
