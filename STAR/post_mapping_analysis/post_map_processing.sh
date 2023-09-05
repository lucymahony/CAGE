#!/bin/bash 

source package 638df626-d658-40aa-80e5-14a275b7464b # source samtools 
source package 6394519c-541f-479c-b064-dd0b912eac04 # bedtools


process_file() {
    tissue=$1
    suffix=$2
    prefix=$3
    intermediate_name=${tissue}_${prefix}
    summary_statistics=${intermediate_name}_summary_statistics.txt

    file=../star_output/${tissue}${suffix}
    echo "Processing: ${file}" > ${summary_statistics}
    echo "---------------------------------------------" >> ${summary_statistics}

    # sort 

    samtools sort --threads 16 ${file} -o ${intermediate_name}_sorted.bam
    # human readable .sam files 
    samtools view -h --threads 16  ${intermediate_name}_sorted.bam > ${intermediate_name}.sam
    echo " Total number of reads in the BAM file " >> ${summary_statistics}
    samtools view -c ${intermediate_name}_sorted.bam >> ${summary_statistics}
    echo " counting only mapped (primary aligned) reads" >> ${summary_statistics}
    samtools view -c -F 260 ${intermediate_name}_sorted.bam >> ${summary_statistics}

    bedtools genomecov -ibam ${intermediate_name}.sam -5 -strand +  -bg > ${intermediate_name}_plus.bed
    bedtools genomecov -ibam ${intermediate_name}.sam -5 -strand -  -bg > ${intermediate_name}_minus.bed 

    echo "Statistics for ${intermediate_name}_plus.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}_plus.bed >> ${summary_statistics}
    echo "Statistics for ${intermediate_name}_minus.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}_minus.bed >> ${summary_statistics}

    cat ${intermediate_name}_plus.bed | awk '{print $1 "\t" $2 "\t" "+" "\t" $4}' > ${intermediate_name}_plus.ctss.bed
    cat ${intermediate_name}_minus.bed | awk '{print $1 "\t" $2 "\t" "-" "\t" $4}' > ${intermediate_name}_minus.ctss.bed
    
    echo "Statistics for  ${intermediate_name}_plus.ctss.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}_plus.ctss.bed >> ${summary_statistics}
    echo "Statistics for  ${intermediate_name}_minus.ctss.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}_minus.ctss.bed >> ${summary_statistics}

    cat ${intermediate_name}_plus.ctss.bed ${intermediate_name}_minus.ctss.bed > ${intermediate_name}.ctss.bed
    echo "Statistics for ${intermediate_name}.ctss.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}.ctss.bed >> ${summary_statistics}

    sort -k1,1 -k2,2n ${intermediate_name}.ctss.bed > ${intermediate_name}_sorted.ctss.bed
    echo "Statistics for ${intermediate_name}_sorted.ctss.bed:" >> ${summary_statistics}
    wc -l ${intermediate_name}_sorted.ctss.bed >> ${summary_statistics}
    echo "---------------------------------------------" >> ${summary_statistics}
    echo "" >> ${summary_statistics}

    # Cleanup
    rm ${intermediate_name}_sorted.bam
    rm ${intermediate_name}_plus.bed
    rm ${intermediate_name}_minus.bed
    rm ${intermediate_name}_plus.ctss.bed
    rm ${intermediate_name}_minus.ctss.bed
    rm ${intermediate_name}.ctss.bed
}

mapfile -t TASKS < tissue.list

# Use the $SLURM_PROCID variable to access the correct tissue and suffix from the list
IFS=' ' read -ra ARGS <<< "${TASKS[$SLURM_PROCID]}"
process_file "${ARGS[0]}" "${ARGS[1]}" "${ARGS[2]}"


