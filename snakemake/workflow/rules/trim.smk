# Trim reads, write number of reads to count file, perform fastqc on trimmed reads 
import pandas as pd

def get_solexa(wildcards):
    """
    Returns the solexa sequencing primers file if the genome is 'fielder'
    """
    sample_name = wildcards.sample.split('_R')[0] # Removes the bit in the case of LE1_R1 e.c.t 
    samples_df = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv', delim_whitespace=True, dtype=str, comment="#", index_col=False)
    # Check if the sample_name exists in the DataFrame
    matching_rows = samples_df.loc[samples_df['sample_name'] == sample_name]
    if matching_rows.empty:
        raise ValueError(f"No entry found for sample {sample_name} in samples.tsv")
    
    genome = samples_df.loc[samples_df['sample_name'] == sample_name, 'genome'].values[0]
    if genome == 'fielder':
        return "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/fielder_cage_data/solexa_sequencing_primers.fa"
    else: 
        return None


rule trim_reads_pe:
    output:
        read_1="{outdir}/{sample}_1P.pe.trim.fastq.gz",
        read_2="{outdir}/{sample}_2P.pe.trim.fastq.gz",
        read_1_unpaired="{outdir}/{sample}_1U.pe.trim.fastq.gz",
        read_2_unpaired="{outdir}/{sample}_2U.pe.trim.fastq.gz"
    input: 
        get_fastqs,
        solexa = get_solexa,
    shell: 
        r"""
        source package 50fcf79b-73a3-4f94-9553-5ed917823423 #  trimmomatic 0.39  
        echo {wildcards.sample} "!!!!" &&
        if [ -n "{input.solexa}" ]; then
            trimmomatic PE -threads {threads} {input[0]} {input[1]} \
            {output.read_1} {output.read_1_unpaired} {output.read_2} {output.read_2_unpaired} \
            ILLUMINACLIP:{input.solexa}:2:30:12 SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20
        else
            trimmomatic PE -threads {threads} {input[0]} {input[1]} \
            {output.read_1} {output.read_1_unpaired} {output.read_2} {output.read_2_unpaired} \
            SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20
        fi
        """


rule trim_reads_se:       
    output: "{outdir}/{sample}.se.trim.fastq.gz"
    input: get_fastqs,
    shell: 
        r"""
        source package 50fcf79b-73a3-4f94-9553-5ed917823423 #  trimmomatic 0.39  
        trimmomatic SE -threads {threads} {input} {output} SLIDINGWINDOW:4:20 MINLEN:20 AVGQUAL:20 
        """


rule trim_count_reads_pe:
    output: temp("{outdir}/{sample}.pe.trim.count")
    input:  
        one_p = "{outdir}/{sample}_1P.pe.trim.fastq.gz",
        one_u = "{outdir}/{sample}_1U.pe.trim.fastq.gz",
        two_p = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
        two_u = "{outdir}/{sample}_2U.pe.trim.fastq.gz",

    shell:
        r"""
        echo {wildcards.sample} 1P: $(( $(wc -l <{input.one_p}) / 4 )) > {output} &&\
        echo {wildcards.sample} 1U: $(( $(wc -l <{input.one_u}) / 4 )) >> {output} &&\
        echo {wildcards.sample} 2P: $(( $(wc -l <{input.two_p}) / 4 )) >> {output} &&\
        echo {wildcards.sample} 2U: $(( $(wc -l <{input.two_u}) / 4 )) >> {output}
        """


rule trim_count_reads_se:
    output: temp("{outdir}/{sample}.se.trim.count")
    input: "{outdir}/{sample}.se.trim.fastq.gz"
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
        echo {wildcards.sample} R1 : $(( $(gunzip {input[0]} | wc -l) / 4 )) > {output} &&\
        echo {wildcards.sample} R2 : $(( $(gunzip {input[1]} | wc -l) / 4 )) >> {output}
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

