log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(ggplot2)
library(reshape2)
library(dplyr)


data <-read.csv(snakemake@input, sept = '\t', header = TRUE)
colnames(data)



# # Melt the data to long format
# data_long <- melt(data, id.vars = "Tissue", variable.name = "Method", value.name = "Value")
# data_long <- data_long %>%
#   mutate(Type = gsub("\\..*", "", Method), 
#          Bar = gsub(".*\\.", "", Method))
# # Create a new variable for the grouped tissue type
# data_long$Tissue_Group <- substr(data_long$Tissue, 1, 2)

# # Calculate means
# data_means <- data_long %>%
#   group_by(Tissue_Group, Type, Bar) %>%
#   summarize(Mean = mean(Value, na.rm = TRUE))

# # Create the bar plot
# p1 <- ggplot(data_means, aes(x = Type, y = Mean, fill = Tissue_Group)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_jitter(data = data_long, aes(x = Type, y = Value, fill = Tissue_Group), 
#                color = "grey", position = position_jitterdodge()) +
#   facet_grid(~Bar) +
#   labs(x = "Method", y = "Value", fill = "Tissue Group") +
#   theme_bw()


# p2 <- ggplot(data_means, aes(x = interaction(Type, Tissue_Group), y = Mean, fill = Bar)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_jitter(data = data_long, aes(x = interaction(Type, Tissue_Group), y = Value), 
#                color = "grey", position = position_jitterdodge()) +
#   scale_fill_manual(values = c("grey", "white"), guide = guide_legend(title = "Type")) +
#   labs(x = "Method and Tissue Group", y = "Value", fill = "Bar") +
#   theme_bw()


# # Create the bar plot
# p3 <- ggplot(data_means, aes(x = interaction(Type, Tissue_Group), y = Mean, fill = Bar)) + 
#   geom_bar(stat = "identity", position = "stack", aes(fill = Tissue_Group)) +
#   geom_jitter(data = data_long, aes(x = interaction(Type, Tissue_Group), y = Value), 
#                color = "grey", position = position_jitterdodge()) +
#   labs(x = "Method and Tissue Group", y = "Value", fill = "Tissue Group") +
#   theme_bw()

# # Create the bar plot
# p4 <- ggplot(data_means, aes(x = interaction(Type, Tissue_Group), y = Mean, fill = Bar)) + 
#   geom_bar(stat = "identity", position = "stack", aes(fill = Tissue_Group), colour = "black") +
#   geom_jitter(data = data_long, aes(x = interaction(Type, Tissue_Group), y = Value), 
#                color = "grey", position = position_jitterdodge()) +
#   labs(x = "Method and Tissue Group", y = "Value", fill = "Tissue Group") +
#   theme_bw()

# # Save bar plot with snakemake@output

# ggsave(filename = snakemake@output, width = 7, height = 7)
