#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH -c 4
#SBATCH -p ei-medium 
#SBATCH -o streme.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=10G # 256

# === CONFIGURATION ===
GENOME="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_assembly.fa"  # <-- update with your genome FASTA
GENOME_SIZES="genome.chrom.sizes"
NUM_REGIONS=21000
REGION_LENGTH=200
OUT_PREFIX="random_regions"
PROMOTERS="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Wheat_expression_prediction/intermediate_data/promoter_seqs.fa"

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0

if [ ! -f "${GENOME}.fai" ]; then
    echo "Indexing genome with samtools faidx..."
    source package c92263ec-95e5-43eb-a527-8f1496d56f1a  # Samtools - 1.18
    samtools faidx "$GENOME"
fi

# === Step 1: Generate genome size file for bedtools ===
cut -f1,2 "${GENOME}.fai" > "$GENOME_SIZES"

# === Step 2: Generate random regions ===
bedtools random -l $REGION_LENGTH -n $NUM_REGIONS -g "$GENOME_SIZES" > "${OUT_PREFIX}.bed"

# === Step 3: Add random strands (+/-) ===
awk 'BEGIN{srand()} {strand = (rand() < 0.5) ? "+" : "-"; print $0"\t.\t"strand}' "${OUT_PREFIX}.bed" > "${OUT_PREFIX}_stranded.bed"

# === Step 4: Extract strand-specific FASTA sequences ===
bedtools getfasta -fi "$GENOME" -bed "${OUT_PREFIX}_stranded.bed" -s -name -fo "${OUT_PREFIX}.fa"

echo "Random sequences written to: ${OUT_PREFIX}.fa"

# === Step 5: Run STREME

singularity exec memesuite.image streme --p ${PROMOTERS} --n ${OUT_PREFIX}.fa     

