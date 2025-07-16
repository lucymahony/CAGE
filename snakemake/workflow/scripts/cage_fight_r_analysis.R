
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
  library(dplyr)
  library(glue)
})

setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data")

setwd("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/")
bw_plus <- BigWigFileList(c("CL1.CS.se.star_plus.unique.bw", "CL2.CS.se.star_plus.unique.bw", "CL3.CS.se.star_plus.unique.bw", "CL4.CS.se.star_plus.unique.bw", "CL5.CS.se.star_plus.unique.bw",
 "CR1.CS.se.star_plus.unique.bw", "CR2.CS.se.star_plus.unique.bw", "CR3.CS.se.star_plus.unique.bw", "CR4.CS.se.star_plus.unique.bw", "CR5.CS.se.star_plus.unique.bw",
 "SIS1.CS.se.star_plus.unique.bw", "SIS2.CS.se.star_plus.unique.bw", "SIS3.CS.se.star_plus.unique.bw",
 "SLE1.CS.se.star_plus.unique.bw", "SLE2.CS.se.star_plus.unique.bw", "SLE3.CS.se.star_plus.unique.bw",
 "SRO1.CS.se.star_plus.unique.bw", "SRO2.CS.se.star_plus.unique.bw", "SRO3.CS.se.star_plus.unique.bw",
 "SSP1.CS.se.star_plus.unique.bw", "SSP2.CS.se.star_plus.unique.bw", "SSP3.CS.se.star_plus.unique.bw"))
bw_minus <- BigWigFileList(c("CL1.CS.se.star_minus.unique.bw", "CL2.CS.se.star_minus.unique.bw", "CL3.CS.se.star_minus.unique.bw", "CL4.CS.se.star_minus.unique.bw", "CL5.CS.se.star_minus.unique.bw",
"CR1.CS.se.star_minus.unique.bw", "CR2.CS.se.star_minus.unique.bw", "CR3.CS.se.star_minus.unique.bw", "CR4.CS.se.star_minus.unique.bw", "CR5.CS.se.star_minus.unique.bw",
"SIS1.CS.se.star_minus.unique.bw", "SIS2.CS.se.star_minus.unique.bw", "SIS3.CS.se.star_minus.unique.bw",
"SLE1.CS.se.star_minus.unique.bw", "SLE2.CS.se.star_minus.unique.bw", "SLE3.CS.se.star_minus.unique.bw",
"SRO1.CS.se.star_minus.unique.bw", "SRO2.CS.se.star_minus.unique.bw", "SRO3.CS.se.star_minus.unique.bw",
"SSP1.CS.se.star_minus.unique.bw", "SSP2.CS.se.star_minus.unique.bw", "SSP3.CS.se.star_minus.unique.bw"))
names(bw_minus) <- c("CL1", "CL2", "CL3", "CL4", "CL5", "CR1", "CR2", "CR3", "CR4", "CR5", "SIS1", "SIS2", "SIS3", "SLE1", "SLE2", "SLE3", "SRO1", "SRO2", "SRO3", "SSP1", "SSP2", "SSP3")
names(bw_plus) <- c("CL1", "CL2", "CL3", "CL4", "CL5" , "CR1", "CR2", "CR3", "CR4", "CR5", "SIS1", "SIS2", "SIS3", "SLE1", "SLE2", "SLE3", "SRO1", "SRO2", "SRO3", "SSP1", "SSP2", "SSP3")

# Extract and print sequence names and lengths for each BigWigFile

print(seqnames(seqinfo(bw_plus)))
print(seqlengths(seqinfo(bw_plus)))
print(seqnames(seqinfo(bw_minus)))
print(seqlengths(seqinfo(bw_minus)))


CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome = seqinfo(bw_plus[[1]])) 
CTSSs <- calcTPM(CTSSs, inputAssay="counts", outputAssay="TPM", outputColumn="subsetTags")

CTSSs <- calcSupport(CTSSs, 
                            inputAssay="counts", 
                            outputColumn="support", 
                            unexpressed=0)
# Remove only expressed in 1 
CTSSs <- subset(CTSSs, support > 1)
CTSSs <- calcTPM(CTSSs, totalTags="totalTags")
CTSSs <- calcPooled(CTSSs)                            
simple_TCs <- clusterUnidirectionally(CTSSs, 
                                     pooledCutoff=0.1, 
                                     mergeDist=20)
                                     
BCs <- clusterBidirectionally(CTSSs, balanceThreshold=0.95)
enhancers <- subset(enhancers, bidirectionality > 0)

supported_enhancers <- subsetBySupport(exampleBidirectional,
                                       inputAssay="counts",
                                       unexpressed=0,
                                       minSamples=1)

simple_TCs <- calcTPM(simple_TCs, 
                                 totalTags = "totalTags")
simple_TCs <- subsetBySupport(simple_TCs,
                                         inputAssay="TPM",
                                         unexpressed=1,minSamples=2)

