# CAGE

This repository contains a CAGE-seq data analysis pipeline.

## Directory Structure

```
CAGE/
├── input_data/
│   ├── *_cage_data/           # Raw CAGE-seq inputs
│   └── *_genome_data/         # Reference genomes and annotations

├── intermediate_data/         

├── results/
│   ├── figs/                  
│   └── data/                 

├── snakemake/
│   └── config/
│       ├── config.yaml       
│       ├── samples.tsv       
│       └── units.tsv         

└── workflow/
    ├── envs/                  
    ├── rules/                
    └── scripts/            
```

## Getting Started

1. Clone the repository:
   ```
   git clone https://github.com/lucymahony/CAGE.git
   ```
