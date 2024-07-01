# Trim reads, write number of reads to count file, perform fastqc on trimmed reads 
include: "common.smk"


rule trim_reads_pe:
    output:
        read_1="{outdir}/{sample}_1P.trim.fastq.gz",
        read_2="{outdir}/{sample}_2P.trim.fastq.gz",
        read_1_unpaired="{outdir}/{sample}_1U.trim.fastq.gz",
        read_2_unpaired="{outdir}/{sample}_2U.trim.fastq.gz"
    input: 
        get_fastqs,
        solexa = "../../../../../tmp/test_CAGE_data/solexa_sequencing_primers.fa",
    conda: "../envs/trim.yaml"
    shell: 
        r"""

        echo {wildcards.sample} "!!!!" &&
        trimmomatic PE -threads {threads} {input[0]} {input[1]} \
        {output.read_1} {output.read_1_unpaired} {output.read_2} {output.read_2_unpaired} \
        ILLUMINACLIP:{input.solexa}:2:30:12 SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20 
        """


rule trim_reads_se:       
    output: "{outdir}/{sample}.trim.fastq.gz"
    input: get_fastqs,
    conda: "../envs/trim.yaml"
    shell: 
        r"""
        trimmomatic SE -threads {threads} {input} {output} SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20 
        """

rule trim_count_reads_pe:
    output: temp("{outdir}/{sample}.pe.trim.count")
    input:  
        one_p = "{outdir}/{sample}_1P.trim.fastq.gz",
        one_u = "{outdir}/{sample}_1U.trim.fastq.gz",
        two_p = "{outdir}/{sample}_2P.trim.fastq.gz",
        two_u = "{outdir}/{sample}_2U.trim.fastq.gz",

    shell:
        r"""
        echo {wildcards.sample} 1P: $(( $(wc -l <{input.one_p}) / 4 )) > {output} &&\
        echo {wildcards.sample} 1U: $(( $(wc -l <{input.one_u}) / 4 )) >> {output} &&\
        echo {wildcards.sample} 2P: $(( $(wc -l <{input.two_p}) / 4 )) >> {output} &&\
        echo {wildcards.sample} 2U: $(( $(wc -l <{input.two_u}) / 4 )) >> {output}
        """


rule trim_count_reads_se:
    output: temp("{outdir}/{sample}.se.trim.count")
    input: "{outdir}/{sample}.trim.fastq.gz"
    shell:
        r"""
        echo {wildcards.sample}: $(( $(wc -l <{input}) / 4 )) > {output}
        """


rule untrim_count_reads_pe:
    output: temp("{outdir}/{sample}.pe.untrim.count"),
    input:
        get_fastqs,
    shell:
        r"""
        echo {wildcards.sample} R1 : $(( $(wc -l <{input[0]}) / 4 )) > {output} &&\
        echo {wildcards.sample} R2 : $(( $(wc -l <{input[1]}) / 4 )) >> {output}
        """


rule untrim_count_reads_se:
    output: temp("{outdir}/{sample}.se.untrim.count")
    input: get_fastqs,
    shell:
        r"""
        echo {wildcards.sample}: $(( $(wc -l <{input}) / 4 )) > {output}
        """


rule all_counts: 
    """
    Concatenates all the count files into a text file. 
    The output of the count_reads rules are temp() and therefore deleted after this rule is enacted.
    """
    output:
        untrimmed = "{outdir}/untrimmed_counts_concatenated.txt",
        trimmed   = "{outdir}/trimmed_counts_concatenated.txt",
    input:
        untrimmed = getAllUntrimmedCountFiles,
        trimmed = getAllTrimmedCountFiles,
    shell:
        r"""
        cat {input.untrimmed} > {output.untrimmed} &&
        cat {input.trimmed} > {output.trimmed}
        """

        

rule fastqc:
    # Runs the fastqc on the now trimmed files 
    # At the moment fastqc is installed in the conda env python_env. Change this at some point!
    # Example command is (python_env) [u10093927@by0q4n] rules $ snakemake ../../tmp/{RO,LE,SP,IS}{1,2}_{1,2}{P,U}.trim.fastqc.html -s trim.smk --use-conda --conda-frontend conda --cores 16      
    input:  "{outdir}/{sample}.trim.fastq.gz"
    output:
        html="{outdir}/{sample}.trim.fastqc.html",
        zipped="{outdir}/{sample}.trim.fastqc.zip",
    shell:
        """
        # Create temporary directory
        tempdir=$(mktemp -d)
        trap "rm -rf $tempdir" EXIT

        fastqc --outdir $tempdir {input}

        # Move outputs to proper location
        html_path="$tempdir/{wildcards.sample}.trim_fastqc.html"
        zip_path="$tempdir/{wildcards.sample}.trim_fastqc.zip"
        mv $html_path {output.html}
        mv $zip_path {output.zipped}
        """

