rule test_rscript_works:
    input: "input_data.csv"
    output: "output_data.csv"
    singularity:
        "cagesoftwareimage.img
    shell:
        "Rscript plot_alignment_results.R {input} {output}"