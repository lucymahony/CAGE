#!/bin/bash 
#SBATCH -p ei-medium
#SBATCH -o speed.out
#SBATCH -c 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mahony@nbi.ac.uk
#SBATCH --mem=256G #56G # 256
#SBATCH --time=2-00:21:00


OUTPUT_DIRECTORY="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/no_merge_reps"


#for SAMPLE in CL_TC.bed SSP_TC.bed SIS_TC.bed SRO_TC.bed SLE_TC.bed;
for SAMPLE in ALL_TC.bed;
do
  singularity exec feb_252.img Rscript speed_up_annotate_bed_tcs.R $OUTPUT_DIRECTORY $SAMPLE
  SAMPLE_LABEL="${SAMPLE%%_*}"  # CR from CR_TC.bed
  CSV_FILE="${OUTPUT_DIRECTORY}/${SAMPLE_LABEL}_annotated.csv"

  OUTPUT_FILE="${OUTPUT_DIRECTORY}/${SAMPLE%.*}_merged_output.csv"
  tail -n +2 "$OUTPUT_DIRECTORY/$SAMPLE" > bed_tmp.txt
  tail -n +2 "$CSV_FILE" | cut -d',' -f5,6 > annot_tmp.txt
  {
    echo "seqnames,start,end,score,strand,dominant_start,dominant_end,annotation,nearest_gene"
    paste bed_tmp.txt annot_tmp.txt |
      awk -F'\t' 'BEGIN {OFS=","} {print $1,$2,$3,$5,$6,$7,$8,$9,$10}'
  } > "$OUTPUT_FILE"
  rm bed_tmp.txt annot_tmp.txt

done
