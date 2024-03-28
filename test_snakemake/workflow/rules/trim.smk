# Trim reads, write number of reads to count file, perform fastqc on trimmed reads 

rule trim_reads:
    output:
        read_1="{outdir}/{sample}_1P.trim.fastq.gz",
        read_2="{outdir}/{sample}_2P.trim.fastq.gz",
        read_1_unpaired="{outdir}/{sample}_1U.trim.fastq.gz",
        read_2_unpaired="{outdir}/{sample}_2U.trim.fastq.gz"
    input: 
        read_1 = "../../../test_CAGE_data/{sample}_R1_test.fastq",
        read_2 = "../../../test_CAGE_data/{sample}_R2_test.fastq",
        solexa = "../../../test_CAGE_data/solexa_sequencing_primers.fa",
    conda: "../envs/trim.yaml"
    shell: 
        r"""
        trimmomatic PE -threads {threads} {input.read_1} {input.read_2} \
        {output.read_1} {output.read_1_unpaired} {output.read_2} {output.read_2_unpaired} \
        ILLUMINACLIP:{input.solexa}:2:30:12 SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20 
        """


rule trim_count_reads:
    output: "{outdir}/{sample}.trim.count"
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



rule untrim_count_reads:
    output: 
        R1 = "{outdir}/{sample}_R1.untrim.count",
        R2 = "{outdir}/{sample}_R2.untrim.count",
    input:
        R1 = "../../../test_CAGE_data/{sample}_R1_test.fastq",
        R2 = "../../../test_CAGE_data/{sample}_R2_test.fastq",
    shell:
        r"""
        echo {wildcards.sample} R1 : $(( $(wc -l <{input.R1}) / 4 )) > {output.R1} &&\
        echo {wildcards.sample} R2 : $(( $(wc -l <{input.R2}) / 4 )) > {output.R2}
        """


TISSUES = ['LE', 'RO', 'SP', 'IS']

rule all_counts:
    # HERE you can see when using expand the outdir has to be in double brackets as this is not an expanded thing. 
    # HERE you can also see how runnning this one rule causes all of the trimming to be performed. 
   output:
       untrimmed = "{outdir}/untrimmed_counts_concatenated.txt",
       trimmed   = "{outdir}/trimmed_counts_concatenated.txt",
   input:
       untrimmed = expand("{{outdir}}/{tissue}{repeats}_R{reads}.untrim.count", tissue=TISSUES, repeats=["1", "2", "3"], reads=["1", "2"]),
       trimmed = expand("{{outdir}}/{tissue}{repeats}.trim.count", tissue=TISSUES, repeats=["1", "2", "3"]),
   shell:
       "cat {input.untrimmed} > {output.untrimmed} ; cat {input.trimmed} > {output.trimmed}"


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

