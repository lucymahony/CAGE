#!/bin/bash 
#SBATCH -p ei-largemem
#SBATCH -o test.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=560G # 256
#SBATCH --time=0-01:42:00




source ~/.bashrc 
mamba activate /hpc-home/mahony/miniforge3
conda activate snakemake

BASE_DIR="../intermediate_data/snakemake_intermediate_data"


# Filter BED files to remove unwanted chromosomes
awk 'NR==FNR {a[$1]; next} $1 in a' $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes $BASE_DIR/SIS1.CS.se.star_plus.unique.bed > $BASE_DIR/SIS1.CS.se.star.plus.unique.bed.filtered
awk 'NR==FNR {a[$1]; next} $1 in a'  $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes $BASE_DIR/SIS1.CS.se.star_minus.unique.bed > $BASE_DIR/SIS1.CS.se.star.minus.unique.bed.filtered

# Check that the filtered BED files are not empty
if [ ! -s $BASE_DIR/SIS1.CS.se.star.plus.unique.bed.filtered ]; then
    echo "Error: plus strand BED file is empty"
    exit 1
fi
if [ ! -s $BASE_DIR/SIS1.CS.se.star.minus.unique.bed.filtered ]; then
    echo "Error: minus strand BED file is empty"
    exit 1
fi
# Print that the filtering is complete
echo "Filtering complete"





# Generate coverage and convert to BigWig
bedtools genomecov -bg -i  $BASE_DIR/SIS1.CS.se.star.plus.unique.bed.filtered -g $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes > $BASE_DIR/SIS1.CS.se.star_plus.unique.bed.bg
sort -k1,1 -k2,2n $BASE_DIR/SIS1.CS.se.star_plus.unique.bed.bg > $BASE_DIR/SIS1.CS.se.star_plus.unique.bed.sorted.bg
bedGraphToBigWig  $BASE_DIR/SIS1.CS.se.star_plus.unique.bed.sorted.bg $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes $BASE_DIR/SIS1.CS.se.star_plus.unique.bw
echo "Plus strand coverage complete"

bedtools genomecov -bg -i  $BASE_DIR/SIS1.CS.se.star.minus.unique.bed.filtered -g $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes > $BASE_DIR/SIS1.CS.se.star_minus.unique.bed.bg
sort -k1,1 -k2,2n $BASE_DIR/SIS1.CS.se.star_minus.unique.bed.bg > $BASE_DIR/SIS1.CS.se.star_minus.unique.bed.sorted.bg
bedGraphToBigWig  $BASE_DIR/SIS1.CS.se.star_minus.unique.bed.sorted.bg $BASE_DIR/Taestivum.ChineseSpring.chrom.sizes $BASE_DIR/SIS1.CS.se.star_minus.unique.bw
echo "Plus strand coverage complete"