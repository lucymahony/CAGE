# CAGE
Scripts in CAGE-seq data analysis pipeline.

Files stored in the following structure:

├── input_data/
│ ├── *_cage_data/
│ ├── *_genome_data/
│
├── intermediate_data/
│
├── results/
│
├── snakemake/
│ ├── config/
│ ││ ├──config.yaml 
│ ││ ├──samples.tsv
│ ││ ├──units.tsv
│ ├── workflow/
│ │├envs/
│ │├rules/
│ │├scripts/














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
