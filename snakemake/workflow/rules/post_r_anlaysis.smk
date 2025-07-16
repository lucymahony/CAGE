rule aggregate_tag_clusters:
    # This script is built off bed_tools_to_create_aggregate.sh
    input:
        expand("{input_dir}/{file}", input_dir="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data", file=[
            "SSP1_TC.bed", "SSP2_TC.bed", "SSP3_TC.bed",
            "SRO1_TC.bed", "SRO2_TC.bed", "SRO3_TC.bed",
            "SLE1_TC.bed", "SLE2_TC.bed", "SLE3_TC.bed",
            "SIS1_TC.bed", "SIS2_TC.bed", "SIS3_TC.bed",
            "CR1_TC.bed", "CR2_TC.bed", "CR3_TC.bed", "CR4_TC.bed", "CR5_TC.bed",
            "CL1_TC.bed", "CL2_TC.bed", "CL3_TC.bed", "CL4_TC.bed", "CL5_TC.bed"
        ])
    output:
        "{input_dir}/shared_{category}_at_least_{number_shared_files}_overlaps.bed",
        "{input_dir}/shared_{category}_at_least_{number_shared_files}_overlaps_max_score.bed"
    params:
        input_dir="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data",
        number_shared_files=6
    shell:
        """
        source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0

        declare -A file_sets
        file_sets["SP_fielder"]="SSP1_TC.bed SSP2_TC.bed SSP3_TC.bed"
        file_sets["RO_fielder"]="SRO1_TC.bed SRO2_TC.bed SRO3_TC.bed"
        file_sets["LE_fielder"]="SLE1_TC.bed SLE2_TC.bed SLE3_TC.bed"
        file_sets["IS_fielder"]="SIS1_TC.bed SIS2_TC.bed SIS3_TC.bed"
        file_sets["R_cadenza"]="CR1_TC.bed CR2_TC.bed CR3_TC.bed CR4_TC.bed CR5_TC.bed"
        file_sets["L_cadenza"]="CL1_TC.bed CL2_TC.bed CL3_TC.bed CL4_TC.bed CL5_TC.bed"
        file_sets['fielder_cadenza_all']="shared_RO_fielder_at_least_2_overlaps_max_score.bed shared_LE_fielder_at_least_2_overlaps_max_score.bed shared_IS_fielder_at_least_2_overlaps_max_score.bed shared_SP_fielder_at_least_2_overlaps_max_score.bed shared_R_cadenza_at_least_2_overlaps_max_score.bed shared_L_cadenza_at_least_2_overlaps_max_score.bed"

        for category in "${!file_sets[@]}"; do
            files=(${file_sets[$category]})
            output_file="{params.input_dir}/shared_${{category}}_at_least_{params.number_shared_files}_overlaps.bed"
            max_score_file="{params.input_dir}/shared_${{category}}_at_least_{params.number_shared_files}_overlaps_max_score.bed"

            final_files=()

            for strand in "+" "-"; do
                temp_files=()
                for file in "${{files[@]}}"; do
                    input_file="{params.input_dir}/${{file}}"
                    temp_file="${{input_file%.bed}}_strand_${{strand}}.bed"
                    awk -v strand="${{strand}}" '$6 == strand' "${{input_file}}" > "${{temp_file}}"
                    temp_files+=("${{temp_file}}")
                done

                combined_file="{params.input_dir}/combined_${{strand}}_${{category}}.bed"
                cat "${{temp_files[@]}}" | sort -k1,1 -k2,2n > "${{combined_file}}"

                if [ ! -s "${{combined_file}}" ]; then
                    echo "No regions found on strand ${{strand}} for category ${{category}}. Skipping..."
                    rm "${{temp_files[@]}}" "${{combined_file}}"
                    continue
                fi

                shared_file="{params.input_dir}/shared_${{category}}_strand_${{strand}}.bed"
                bedtools intersect -a "${{combined_file}}" -b "${{temp_files[@]}}" -c | \
                    awk -v threshold="{params.number_shared_files}" '$NF >= threshold' > "${{shared_file}}"

                if [ ! -s "${{shared_file}}" ]; then
                    echo "No shared regions found on strand ${{strand}} for category ${{category}}. Skipping..."
                    rm "${{temp_files[@]}}" "${{combined_file}}" "${{shared_file}}"
                    continue
                fi

                final_file="{params.input_dir}/final_${{category}}_strand_${{strand}}.bed"
                bedtools merge -i "${{shared_file}}" -d 10 -c 5,6,7,8 -o collapse,distinct,collapse,collapse | \
                    awk 'BEGIN{{OFS="\t"}} {print $1, $2, $3, ".", $4, $5, $6, $7, $8}' > "${{final_file}}"
                final_files+=("${{final_file}}")

                rm "${{temp_files[@]}}" "${{combined_file}}" "${{shared_file}}"
            done

            if [ ${{#final_files[@]}} -gt 0 ]; then
                cat "${{final_files[@]}}" | sort -k1,1 -k2,2n > "${{output_file}}"
                echo "Final shared regions for ${{category}} saved to ${output_file}"
            else
                echo "No shared regions found across all files for category ${{category}}."
                continue
            fi

            rm -f "${{final_files[@]}}"

            awk 'BEGIN{{OFS="\t"}} {{
                split($5, scores, ",");
                split($7, starts, ",");
                split($8, ends, ",");

                max_index = 1;
                max_score = scores[1];
                for (i = 2; i <= length(scores); i++) {{
                    if (scores[i] > max_score) {{
                        max_score = scores[i];
                        max_index = i;
                    }}
                }}

                print $1, $2, $3, $4, max_score, $6, starts[max_index], ends[max_index];
            }}' "${{output_file}}" > "${{max_score_file}}"

            echo "Final shared regions with max score for ${{category}} saved to ${max_score_file}"
        done
        """