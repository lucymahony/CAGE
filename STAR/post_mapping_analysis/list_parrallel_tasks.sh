#!/bin/bash 

declare -a tissues=("IS1" "IS2" "IS3" "SP1" "SP2" "SP3" "RO1" "RO2" "RO3" "LE1" "LE2" "LE3") 
declare -a input_sam_name=("uniquelymapped_${tissue}.sam" "${tissue}.out.sam")
declare -a prefixes=("UO" "UM")

> tissue.list
for tissue in "${tissues[@]}"; do
    for index in "${!input_sam_name[@]}"; do
        echo "${tissue} ${input_sam_name[$index]} ${prefixes[$index]}" >> tissue.list
   done
done


