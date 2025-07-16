suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
input_file_path <- args[1]
output_file_path <- args[2]

process_summary_statistics_text_file <- function(input_file_path) {
  text <- readLines(input_file_path)
  total_reads <- c()
  primary_mapped <- c()
  sample_names <- c()
  current_total <- NA
  current_primary <- NA
  for (i in seq_along(text)) {
    line <- text[i]
    if (grepl("Total number of reads", line)) {
      current_total <- as.numeric(text[i + 1])
    }
    
    if (grepl("Counting only mapped", line)) {
      current_primary <- as.numeric(text[i + 1])
    }
    
    if (grepl("Statistics for", line)) {
      sample_name <- gsub(".*intermediate_data/(.*)\\..*", "\\1", line)
      sample_name <- gsub("\\..*", "", sample_name)
      sample_names <- c(sample_names, sample_name)
      total_reads <- c(total_reads, current_total)
      primary_mapped <- c(primary_mapped, current_primary)
    }
  }
  df <- data.frame(
    Sample_Name = sample_names,
    Total_Reads = total_reads,
    Primary_Mapped_Reads = primary_mapped
  )
  df <- df[!duplicated(df), ]
  return(df)
}

create_mapping_plot <- function(df) {

  # Add a tissue type column based on Sample_Name
  df <- df %>%
    mutate(
      Tissue = case_when(
        grepl("^IS", Sample_Name) ~ "Fielder Immature Spike PE",
        grepl("^SP", Sample_Name) ~ "Fielder Spike PE",
        grepl("^RO", Sample_Name) ~ "Fielder Root PE",
        grepl("^LE", Sample_Name) ~ "Fielder Leaf PE",
        grepl("^SIS", Sample_Name) ~ "Fielder Immature Spike",
        grepl("^SSP", Sample_Name) ~ "Fielder Spike",
        grepl("^SRO", Sample_Name) ~ "Fielder Root",
        grepl("^SLE", Sample_Name) ~ "Fielder Leaf",
        grepl("^CL", Sample_Name) ~ "Cadenza Leaf",
        grepl("^CR", Sample_Name) ~ "Cadenza Root",
        grepl("^CSS", Sample_Name) ~ "Chinese Spring Spike",
        grepl("^CSEM", Sample_Name) ~ "Chinese Spring Embryo",
        grepl("^CSRO", Sample_Name) ~ "Chinese Spring Root",
        TRUE ~ "Unknown"
      )
    )
  df <- df %>%
    filter(!Tissue %in% c("Fielder Immature Spike PE", "Fielder Spike PE", "Fielder Root PE", "Fielder Leaf PE"))
  
  df <- df %>%
    mutate(
      Percent_Primary_Mapped = (Primary_Mapped_Reads / Total_Reads) * 100
    )
  plot <- ggplot(df, aes(x = Tissue, y = Percent_Primary_Mapped, fill = Tissue)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_point(position = position_dodge(width = 0.9), size = 3, alpha=0.9, color='grey') +
    theme_minimal() +
    labs(
      x = "Tissue",
      y = "Percentage of Reads Mapped",
      title = "Percentage of Primary Mapped Reads per Tissue"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(plot)
}

df <- process_summary_statistics_text_file(input_file_path)
print(df)
plot <- create_mapping_plot(df)
ggsave(filename = output_file_path, plot, width = 7, height = 7)
print('Plot saved at ')
