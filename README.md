# CAGE

This repository contains a CAGE-seq data analysis pipeline.


![Snakemake Workflow Diagram](docs/images/snakemake_pipeline.drawio.pdf)
![Snakemake Workflow Diagram](docs/images/Untitled\ Diagram.drawio.png)
## Directory Structure

```
CAGE/
├── input_data/
│   ├── *_cage_data/           # Raw CAGE-seq inputs
│   └── *_genome_data/         # Reference genomes and annotations

├── tmp/         

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
2. File paths:
    Set up file paths in the config file units.tsv copying the file format
    Specify tissue and genome in samples.tsv. Use the genome naming convention that matches the directories in input_data/

3. Running analysis
    Specify desired output files in `snakemake/run_snakemake.sh`
