#!/bin/sh
#SBATCH -p ei-medium
#SBATCH -c 1	
#SBATCH --mem 10G				# memory pool for all cores
#SBATCH --time=0-01:00:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

source package c92263ec-95e5-43eb-a527-8f1496d56f1a  # Samtools - 1.18
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64  # BedTools2 2.30.0

make_small_files () {
    # Define variables
    chr=$1
    start=$2
    stop=$3

    # Input file paths
    reference_assembly=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna 
    reference_annotatation=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3
    
    intermediate_data_directory=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data
    all_cage_tag_clusters=${intermediate_data_directory}/no_merge_reps/ALL_TC.bed
    cadenza_leaf_cage_tag_clusters=${intermediate_data_directory}/no_merge_reps/CL_TC.bed
    
    region_directory="${chr}_${start}_${stop}"
    # Check if directory exists, if not create it - specific to region. 
    if [ ! -d "${intermediate_data_directory}/${region_directory}" ]; then
      mkdir -p "${intermediate_data_directory}/${region_directory}"
    fi

    # Genome assembly - sed '$d' deletes the last line
    sed -n "/^>Chr${chr}\b/,/^>/p" "$reference_assembly" | sed '$d' > "${intermediate_data_directory}/${region_directory}/subsection_IWGSC_CS_RefSeq_v2.1.fna"
    echo "Genome assembly extraction complete"

    # Genome annotation - Subsect using bedtools 
    region_gff="${intermediate_data_directory}/${region_directory}/subsection_annotation_200916_HC.gff3"
    echo -e "Chr${chr}\t${start}\t${stop}" > region.bed
    bedtools intersect -wa -a "$reference_annotatation" -b  region.bed > "$region_gff"
    sed -i '1s/^chr/Chr/' "${region_gff}" # Change chr1A to Chr1A e.g. capitalise first letter of first row
    echo "Genome annotation extraction complete. Output:"

    # Filter CAGE tag clusters for region of interest - Subsect using bedtools
    # ALL 
    
    bedtools intersect -wa -a  "$all_cage_tag_clusters"  -b  region.bed  > "${intermediate_data_directory}/${region_directory}/ALL_CAGE_TCs_my_annot.bed"
    awk '$5 > 0.5' "${intermediate_data_directory}/${region_directory}/ALL_CAGE_TCs_my_annot.bed" > "${intermediate_data_directory}/${region_directory}/FILTERED_ALL_CAGE_TCs_my_annot.bed"
    echo "CAGE tag clusters extraction complete. Output:"
    # CAGE tag clusters extraction complete. Output:
    # filter by 5th column being more than 0.5
    head "${intermediate_data_directory}/${region_directory}/FILTERED_ALL_CAGE_TCs_my_annot.bed"


    ## CL - cadenza_leaf_cage_tag_clusters
    #echo "head of cadenza_leaf_cage_tag_clusters"
    #head "$cadenza_leaf_cage_tag_clusters"
    #bedtools intersect -wa -a  "$cadenza_leaf_cage_tag_clusters"  -b  region.bed  > "${intermediate_data_directory}/${region_directory}/CL_CAGE_TCs_my_annot.bed"
    #awk '$5 > 0.5' "${intermediate_data_directory}/${region_directory}/CL_CAGE_TCs_my_annot.bed" > "${intermediate_data_directory}/${region_directory}/CL_CAGE_TCs_my_annot.bed"
    #echo "CAGE tag clusters extraction complete. Output:"
    #head "${intermediate_data_directory}/${region_directory}/CL_CAGE_TCs_my_annot.bed"
    rm region.bed
    
    # Filter CAGE BAM files for region of interest
    
    #for SAMPLE in SIS1 SIS2 SIS3 SLE1 SLE2 SLE3 SRO1 SRO2 SRO3 SSP1 SSP2 SSP3 CL1 CL2 CL3 CL4 CL5 CR1 CR2 CR3 CR4 CR5
    for SAMPLE in CL1
    do
        cage_bam_file=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/${SAMPLE}.CS.se.star.bam
        cage_section_bam_file=${intermediate_data_directory}/${region_directory}/${SAMPLE}.bam
        cage_sam_file=${intermediate_data_directory}/${region_directory}/${SAMPLE}.sam
        # index and sort BAM file    
        samtools index -c $cage_bam_file 
        samtools view -b "$cage_bam_file" "Chr$chr:$start-$stop" | samtools sort -o "$cage_section_bam_file"
        ## convert to sam file - this can be used in IGV
        samtools view -h $cage_section_bam_file > $cage_sam_file
        rm $cage_section_bam_file
        echo "CAGE SAM file extraction complete"
    done
}

#make_small_files '1A' '1293900' '12941173'
#make_small_files '1B' '1489800' '14900216'
#make_small_files '1D' '1013000' '10145801'

#make_small_files '1A' '13046000' '13050288'
#make_small_files '1B' '15577000' '15581389'
#make_small_files '1D' '10270000' '10274452'
make_small_files '7A'	'738262315' '738267073'
make_small_files '7B'	'754110234' '754112414'
make_small_files '7D'	'639808726' '639814470'
