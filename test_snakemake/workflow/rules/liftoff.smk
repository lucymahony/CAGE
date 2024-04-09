# This is used for generating fielder and Cadenza annotations from CS 
rule liftoff:
    input:
        closely_related_assembly='',
        reference_annotation='',
        reference_assembly='',
        chroms='',
    output: "{outdir}/output.gff"
    shell: r"liftoff {input.closely_related_assembly} {input.reference_assembly} -g {input.reference_annotation} -chroms {input.chroms} -o {output}"
