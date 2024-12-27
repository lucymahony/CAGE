rule rscript_plot_alignment_results:
    input: "{outdir}/{genome}.{mappingtool}.all_post_align_summary_statistics.txt"
    output: "{outdir}/{genome}.{mappingtool}.alignment_statistics_results.pdf"
    singularity:
        "workflow/scripts/cagesoftwareimage.img"
    shell:
        "Rscript workflow/scripts/plot_alignment_results.R {input} {output}"


rule rscript_cage_r_analysis:
    input: 
        CTSSNBedfiles=getAllSortedCTSSNBedfiles,
        input_files_directory="{outdir}", 
        samples_tsv_file_path='/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv', 
        
    params:
        file_pattern="sorted.ctss.n.bed",

    output: "{outdir}/{genome_aligned_to}.{mapping_tool}.cage_r_analysis_results.txt"

    singularity:
        "workflow/scripts/nov_24.img"

    shell:
        "Rscript workflow/scripts/cage_r_analysis.R {params.file_pattern} {input.input_files_directory} {input.samples_tsv_file_path}"





