

library(GenomicRanges)
library(rtracklayer)




gff_file <- "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3"  
print(gff_file)
gff <- rtracklayer::import(gff_file)
print(gff)
gff$transcript_type <- gff$biotype
gff$gene_name <- gff$Name
print(gff)

