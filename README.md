# CAGE
Scripts in CAGEseq data analysis pipeline.


## STAR
Mapping reads using STAR.
### Indexing 
`cd indexing`
Requires the genome sequence and annotation in the folder input_data/fielder_genome_data/ and the empty folder `genome_directory` to store the output. 

1. `sbatch fielder_indexing.sh`

### Mapping 
`cd mapping`

1. `sbatch `

`cd post_mapping_analysis`

1. `sbatch `
