# CAGE

This repository contains a CAGE-seq data analysis pipeline.

## Directory Structure

The project is organized into the following structure:

CAGE/
├── input_data/
│ ├── *_cage_data/ # Raw CAGE-seq inputs
│ └── *_genome_data/ # Reference genomes and annotations
│
├── intermediate_data/ # Processed files and intermediate outputs
│
├── results/
│ ├── figs/ # Figures for plots and visualisations
│ └── data/ # Final result files (e.g., tables)
│
├── snakemake/
│ └── config/
│ ├── config.yaml # Snakemake configuration
│ ├── samples.tsv # Sample metadata
│ └── units.tsv # Unit-level metadata
│
└── workflow/
├── envs/ # Conda environments for each rule
├── rules/ # Snakemake rule definitions
└── scripts/ # Custom Python/R scripts used in rules

## Getting Started

To run the pipeline:

1. Clone the repository:
   ```bash
   git clone https://github.com/lucymahony/CAGE.git
   cd CAGE
