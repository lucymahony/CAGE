### Aligning the reads, first create the tools index files and then align them 
# Possibly the assemby and annotation names should be saved as tuples in the config file?
# For fielder the assemly name is  ____ for fielder the annotation is ________- note I need to add the liftover into the pipeline?
# Index files are called assembly_name.indexer_tool_name.index (maybe using the annotation name would be equivalent?)

rule star_genome_index:
    input:
        fasta = "../../../test_CAGE_data/{assembly_name}",  
        gtf = "../../../test_CAGE_data/{annotation_name}"  
    output: directory("{outdir}/{assembly_name}.star.index")
    shell:
        """
        mkdir -p {output} ;
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 \
            --limitGenomeGenerateRAM 100111047092 \
            --genomeDir {output}        
        """


rule bowtie2_genome_index:
    params:
        basename="{outdir}/bowtie_index_{genome}"
    input: "../../../test_CAGE_data/{assembly_name}"
    output: 
        multiext(
            "{outdir}/{assembly_name}.bowtie2.index.",
            "1.bt2l",
            "2.bt2l",
            "3.bt2l",
            "4.bt2l",
            "rev.1.bt2l",
            "rev.2.bt2l",
         )
    shell: 
        r"""
        bowtie2-build {input} {params.basename} --threads={threads}
        """


rule bwa_genome_index:
    output: 
        multiext(
            "{outdir}/{assembly_name}.bwa.index",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa"), 
    input: "../../../test_CAGE_data/{assembly_name}"
    shell:
        r"""
        bwa index {input} -p {wildcards.outdir}/bwa_index_{wildcards.genome_stem}
        """

rule hisat2_genome_index:
    input:
        fasta = "../../../test_CAGE_data/{assembly_name}",  
        gtf = "../../../test_CAGE_data/{annotation_name}"  

    output: "{outdir}/{assembly_name}.hisat2.index"
    shell:
        r""" 
        hisat2_extract_splice_sites.py {input.gtf} > {wildcards.outdir}/fielder.ss
        hisat2_extract_exons.py {input.gtf} > {wildcards.outdir}/fielder.exon
        hisat2-build -p {threads} --ss {wildcards.outdir}/fielder.ss --exon {wildcards.outdir}/fielder.exon {input.fasta} {output}
        """

