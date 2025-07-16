# Plot the number of mapped reads between Cadenza and CS for each of the cadenza datasets
# Also plot the number of mapped reads between Cadenza and CS for each of the chinese spring datasets

#grep "../intermediate_data/snakemake_intermediate_data/CR4.cadenza.se.star.unique.sorted.ctss.n.bed" cadenza.star.all_post_align_summary_statistics.txt
#Statistics for ../intermediate_data/snakemake_intermediate_data/CR4.cadenza.se.star.unique.sorted.ctss.n.bed:
#1980667 ../intermediate_data/snakemake_intermediate_data/CR4.cadenza.se.star.unique.sorted.ctss.n.bed
#[mahony@EI-HPC interactive snakemake_intermediate_data]$ grep "../intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed" CS.star.all_post_align_summary_statistics.txt
#Statistics for ../intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed:
#1896915 ../intermediate_data/snakemake_intermediate_data/CR4.CS.se.star.unique.sorted.ctss.n.bed

library(ggplot2)
library(reshape2)

#CL1.cadenza  1511339 
#CL2.cadenza  1558728     
#CL3.cadenza  1214893     
#CL4.cadenza  1211526     
#CL5.cadenza  1316006     
#CR1.cadenza  1798300     
#CR2.cadenza  1311285     
#CR3.cadenza  1863797     
#CR4.cadenza  1980667     
#CR5.cadenza  2305261   

#CL1.CS 1438641
#CL2.CS 1495230
#CL3.CS 1153914
#CL4.CS 1148918
#CL5.CS 1247982 
#CR1.CS 1724074
#CR2.CS 1246539
#CR3.CS 1785720
#CR4.CS 1896915
#CR5.CS 2210283
# Data
names <- c("CL1", "CL2", "CL3", "CL4", "CL5", "CR1", "CR2", "CR3", "CR4", "CR5", "CL1", "CL2", "CL3", "CL4", "CL5", "CR1", "CR2", "CR3", "CR4", "CR5")
tissue <- c("Leaf", "Leaf", "Leaf", "Leaf", "Leaf", "Root", "Root", "Root", "Root", "Root", "Leaf", "Leaf", "Leaf", "Leaf", "Leaf", "Root", "Root", "Root", "Root", "Root")
groups <- c("Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Cadenza", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring", "Chinese_Spring")
values <- c(1511339, 1558728, 1214893, 1211526, 1316006, 1798300, 1311285, 1863797, 1980667, 2305261, 1438641, 1495230, 1153914, 1148918, 1247982, 1724074, 1246539, 1785720, 1896915, 2210283)
data <- data.frame(names, groups, tissue, values)

p <- ggplot(data, aes(x=groups, y=values, fill=groups)) +
  geom_point(stat="identity", aes(color=tissue)) +
  geom_line(aes(group=names)) + # add line between matching names
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Reference genome') +
  ylab('Number of uniquley mapped reads') +
  scale_fill_manual(values=c("#009E73", "#D55E00")) + 
  theme_bw()

ggsave("cadenza_vs_CS.pdf", plot=p, width=10, height=10, units="cm", dpi=900)

print("Cadenza")
print(mean(data$values[data$groups == "Cadenza"]))
print("CS")
print(mean(data$values[data$groups == "Chinese_Spring"]))
