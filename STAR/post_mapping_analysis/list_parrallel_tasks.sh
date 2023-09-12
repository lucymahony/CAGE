#!/bin/bash 

declare -a tissues=("IS1" "IS2" "IS3" "SP1" "SP2" "SP3" "RO1" "RO2" "RO3" "LE1" "LE2" "LE3") 
declare -a suffixes=("uniquelymapped_${tissue}.sam" "${tissue}.out.sam")
declare -a prefixes=("UO" "UM")

> tissue.list
for tissue in "${tissues[@]}"; do
    for index in "${!suffixes[@]}"; do
        echo "${tissue} ${suffixes[$index]} ${prefixes[$index]}" >> tissue.list
   done
done


