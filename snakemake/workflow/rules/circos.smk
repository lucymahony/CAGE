import pandas as pd
samples = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)
genome_groups = samples.groupby("genome")["sample_name"].apply(list).to_dict() # A dict where the key is genome name e.g. fielder and the value is a list of the sample names e.g. LE 
genomes = set(samples["genome"].to_list()) # fielder, cadenza, chinese_spring

rule create_genome_bins:
    input: "{outdir}/genome.txt"
    output: "{outdir}/genome_bins.bed"
    shell:
        """
        bedtools makewindows -g {input} -w 1000000 > {output}
        """


rule calculate_average_density:
    input:
        bins="{outdir}/genome_bins.bed",
        beds=lambda wildcards: [f"{wildcards.outdir}/{sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.sorted.ctss.n.bed" for sample in genome_groups.get(wildcards.genome, [])]    
    output:"{outdir}/{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.average_density.bed" # e.g. fielder_mappedto_CS.se.star.average_density.bed
    params:
        combined="{genome_aligned_to}_combined_density.bed"
    shell:
        """
        paste {input.beds} > {params.combined}
        
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
        }}' {params.combined} > {output}
        
        # Clean up intermediate files
        rm {params.combined}
        """



rule generate_circos_plot:
    input:
        karyotype="{outdir}/karyotype.txt",
        circos_conf="{outdir}/circos.conf",
        density_bed_plots=lambda wildcards: [
            f"{wildcards.outdir}/{genome}_mappedto_{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.average_density.bed"
            for genome in genomes]
            
            #if os.path.exists(f"{wildcards.outdir}/{genome}_mappedto_{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.average_density.bed")]
            
    output:"{outdir}/circos_for_{genome}_mappedto_{genome_aligned_to}.{read_type}.{mapping_tool}.png"
    shell:
        """
        echo {input.density_bed_plots}
        source package 22619b3e-43a5-4546-ab14-4561f701f247
        circos -conf {input.circos_conf} 
        """

rule all:
    input:
        "../intermediate_data/snakemake_intermediate_data/circos.png"


#genome.txt
#Chr1A   598660471   
#Chr1B   700547350
#Chr1D   498638509
#Chr2A   787782082
#Chr2B   812755788
#Chr2D   656544405
#Chr3A   754128162 
#Chr3B   851934019
#Chr3D   619618552
#Chr4A   754227511 
#Chr4B   673810255
#Chr4D   518332611
#Chr5A   713360525
#Chr5B   714697677
#Chr5D   569951140  
#Chr6A   622669697
#Chr6B   731188232
#Chr6D   495380293 
#Chr7A   744491536
#Chr7B   764072961
#Chr7D   642921167
# circos.conf
#karyotype = karyotype.txt
#
#<ideogram>
#    <spacing>
#        default = 0.005r
#    </spacing>
#
#    radius           = 0.65r
#    thickness        = 20p
#    fill             = yes
#    stroke_color     = black
#    stroke_thickness = 1p
#
#    show_label       = yes
#    label_font       = bold
#    label_radius     = dims(ideogram,radius_outer) + 0.02r
#    label_size       = 30p
#    label_parallel   = yes
#</ideogram>
#
#<ticks>
#    show_ticks        = yes
#    show_tick_labels  = yes
#
#    <tick>
#        spacing        = 10u
#        size           = 10p
#        thickness      = 1p
#        color          = black
#        show_label     = yes
#        label_size     = 10p
#        label_offset   = 2p
#        format         = %d
#    </tick>
#</ticks>
#
#<plots>
#    # Histogram for density plots
#    <plot>
#        type       = histogram
#        file       = fielder_average_density.bed
#        r0         = 0.7r
#        r1         = 0.9r
#        min        = 0
#        max        = 15
#        color      = blue
#        thickness  = 1
#        fill_under = yes
#        fill_color = blue_a3
#    </plot>
#
#    <plot>
#        type       = histogram
#        file       = cadenza_average_density.bed
#        r0         = 0.5r
#        r1         = 0.7r
#        min        = 0
#        max        = 15
#        color      = red
#        thickness  = 1
#        fill_under = yes
#        fill_color = red_a3
#    </plot>
#
#    <plot>
#        type       = histogram
#        file       = chinese_spring_average_density.bed
#        r0         = 0.3r
#        r1         = 0.5r
#        min        = 0
#        max        = 15
#        color      = green
#        thickness  = 1
#        fill_under = yes
#        fill_color = green_a3
#    </plot>
#</plots>
#
#
#<image>
#    dir           = .
#    file          = circos.png
#    radius        = 1500p
#    background    = white
#    angle_offset  = -90
#    auto_alpha_colors = yes
#    auto_alpha_steps  = 5
#</image>
#
#<<include etc/colors_fonts_patterns.conf>>
#<<include etc/housekeeping.conf>>
#
#