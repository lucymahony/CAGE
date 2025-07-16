import pandas as pd
samples = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)
genome_groups = samples.groupby("genome")["sample_name"].apply(list).to_dict() # A dict where the key is genome name e.g. fielder and the value is a list of the sample names e.g. LE 
genomes = set(samples["genome"].to_list()) # fielder, cadenza, chinese_spring
units = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/units.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)


rule create_genome_bins:
    input: "{outdir}/genome.txt"
    output: "{outdir}/genome_bins.bed"
    shell:
        """
        bedtools makewindows -g {input} -w 1000000 > {output}
        """


def get_bed_files(wildcards):
    samples_given_genome = genome_groups.get(wildcards.genome, [])
    #>>> genome_groups.get('fielder', [])
        #['LE1', 'LE2', 'LE3', 'IS1', 'IS2', 'IS3', 'SP1', 'SP2', 'SP3', 'RO1', 'RO2', 'RO3']
    print(samples_given_genome)
    filtered_samples = []
    if wildcards.read_type == 'se':
        # Filter the samples where the 'number_reads' is 1
        filtered_samples = units[(units['sample'].isin(samples_given_genome)) & (units['number_reads'] == '1')]
    elif wildcards.read_type == 'pe':
        # Filter the samples where the 'number_reads' is 2
        filtered_samples = units[(units['sample'].isin(samples_given_genome)) & (units['number_reads'] == '2')]
    filtered_sample_names = filtered_samples['sample'].tolist()
    return filtered_sample_names


rule calculate_average_density:
    input:
        bins="{outdir}/genome_bins.bed",
        beds=lambda wildcards: [f"{wildcards.outdir}/{sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.unique.sorted.ctss.h.bed" for sample in get_bed_files(wildcards)]    
    output:"{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density.bed" # e.g. fielder_mappedto_CS.se.star.average_density.bed 
    params:
        combined="{genome_aligned_to}_combined_density.bed"
    shell:
        """
        
        echo "Input files:"
        echo {input.beds}
        paste {input.beds} > {wildcards.outdir}/{wildcards.genome}_combined_density.bed
        echo  {wildcards.outdir}/{wildcards.genome}_combined_density.bed

        # Calculate average density
        awk '{{ 
            if (NR == 1) {{
                print "chrom\\tstart\\tend\\taverage_density"
            }} else {{
                sum = 0
                for (i=4; i<=NF; i+=4) sum += $i
                avg = sum / ((NF-3)/4)
                print $1"\\t"$2"\\t"$3"\\t"avg
            }}
        }}' {wildcards.outdir}/{wildcards.genome}_combined_density.bed > {output}
        
        """

rule add_end_average_density:
    input: "{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density.bed"
    output: "{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density_rounded.bed"
    shell:
        """
        awk 'NR==1 {{print; next}} {{$3 = $2 + 1; printf "%s\\t%s\\t%s\\t%.0f\\n", $1, $2, $3, $4}}' {input} > {output}
        
        # Rounds the fourth column ($4) to the nearest integer for circos compatibility
        #awk 'NR==1 {{print; next}} {{printf "%s\\t%s\\t%s\\t%.0f\\n", $1, $2, $3, $4}}' {input} > {output}
        """

rule aggregate_average_density:
    input:
        "{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density_rounded.bed"
    output:
        "{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density_aggregate.bed"
    shell:
        """
        awk '
        BEGIN {{
            BIN_SIZE = 10000000;
        }}
        NR == 1 {{
            print; next;
        }}
        {{
            bin_start = int($2 / BIN_SIZE) * BIN_SIZE;
            bin_end = bin_start + BIN_SIZE;
            bin_key = $1 "\\t" bin_start "\\t" bin_end;
            density_sum[bin_key] += $4;
            count[bin_key] += 1;
        }}
        END {{
            for (key in density_sum) {{
                print key "\\t" density_sum[key] / count[key];
            }}
        }}
        ' {input} > {output}
        """

        

rule generate_circos_plot:
    input:
        karyotype="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/workflow/scripts/karyotype.txt",
        circos_conf="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/workflow/scripts/circos.conf",
        density_bed_plots=lambda wildcards: [
            f"{wildcards.outdir}/{genome}_mappedto_{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.average_density_aggregate.bed"
            for genome in genomes] # genomes = fielder, cadenza, chinese_spring
    output:"{outdir}/circos.{genome_aligned_to}.{read_type}.{mapping_tool}.png"
    shell:
        """
        echo {input.density_bed_plots}
        source package 22619b3e-43a5-4546-ab14-4561f701f247
        circos -conf {input.circos_conf} 
        mv circos.png {output}
        echo "Note the colour platte is Fielder: Blue, Cadenza: Red, Chinese Spring: Green"
        """

