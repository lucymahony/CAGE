### Aligning the reads, first create the tools index files and then align them 
# Possibly the assemby and annotation names should be saved as tuples in the config file?
# For fielder the assemly name is  ____ for fielder the annotation is ________- note I need to add the liftover into the pipeline?
# Index files are called assembly_name.indexer_tool_name.index (maybe using the annotation name would be equivalent?)
# The test files are a small bacterial genome, just used in debugging test files. 


# At some point move these helper functions to the config file?
def get_fasta(variety):
    variety_assembly_dict = {'CS': '../../../../../tmp/test_datasets/iwgsc_refseqv2.1_assembly.fa', 
                    'fielder': '../../../../../tmp/test_datasets/fielder.fa', 
                    'test': '../../../../../tmp/test_datasets/test.fa'}
    return variety_assembly_dict[variety]

def get_gtf(variety):
    variety_annotation_dict = {'CS': '../../../../../tmp/test_datasets/iwgsc_refseqv2.1_annotation_200916_HC.gff3', 
                    'fielder': '../../../../../tmp/test_datasets/fielder.gtf',
                    'test': '../../../../../tmp/test_datasets/test.gff3'}
    return variety_annotation_dict[variety]


########## Genome Indexing ##########

rule star_genome_index:
    # qig snakemake -F ../../../../../tmp/CS.star.index -s align.smk --use-conda --conda-frontend conda --cores 16
    input:
        fasta = lambda wildcards : get_fasta(wildcards.variety),  
        gtf = lambda wildcards : get_gtf(wildcards.variety),
    output: directory("{outdir}/{variety}.star.index")
    shell:
        """
        # Create temporary directory
        tempdir=$(mktemp -d)
        trap "rm -rf $tempdir" EXIT
        
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 \
            --limitGenomeGenerateRAM 100111047092 \
            --genomeDir $tempdir \
        
        # Move outputs to proper location
        mv $tempdir {output} 
        """


rule bowtie2_genome_index:
    input: lambda wildcards : get_fasta(wildcards.variety)
    output: 
        multiext(
            "{outdir}/{variety}.bowtie2.index.",
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",
         )
    shell: 
        r"""
        bowtie2-build {input} {wildcards.outdir}/{wildcards.variety}.bowtie2.index
        """


rule bowtie2_genome_index_large:
    input: lambda wildcards : get_fasta(wildcards.variety)
    output: 
        multiext(
            "{outdir}/{variety}.bowtie2.index.",
            "1.bt2l",
            "2.bt2l",
            "3.bt2l",
            "4.bt2l",
            "rev.1.bt2l",
            "rev.2.bt2l",
         )
    shell: 
        r"""
        bowtie2-build {input} {wildcards.outdir}/{wildcards.variety}.bowtie2.index
        """


rule bwa_genome_index:
    input: lambda wildcards : get_fasta(wildcards.variety)
    output: 
        multiext(
            "{outdir}/{variety}.bwa.index.",
            "amb",
            "ann",
            "bwt",
            "pac",
            "sa"), 
    shell:
        r"""
        bwa index {input} -p {wildcards.outdir}/{wildcards.variety}.bwa.index
        """


rule hisat2_genome_index:
    input:
        fasta = lambda wildcards : get_fasta(wildcards.variety),  
        gtf = lambda wildcards : get_gtf(wildcards.variety),
    output: expand("{{outdir}}/{{variety}}.hisat2.index.{num}.ht2", num=[1, 2, 4, 5, 6, 7, 8])
    shell:
        r""" 
        hisat2_extract_splice_sites.py {input.gtf} > {wildcards.outdir}/fielder.ss
        hisat2_extract_exons.py {input.gtf} > {wildcards.outdir}/fielder.exon
        hisat2-build --ss {wildcards.outdir}/fielder.ss --exon {wildcards.outdir}/fielder.exon {input.fasta} {output}
        """


########## Aligning Paired Trimmed Reads  ##########

rule star_align:
    # snakemake ../../../../../tmp/LE1.test.star.Aligned.sortedByCoord.out.bam -s align.smk --use-conda --conda-frontend conda --cores 16
    # For somereason qig at the start of this line seems to make the conda env not work correctly ??
    # At some point I could come and clean up the output file name to remove the .Aligned.sortedByCoord.out proportion
    input:
        star_index_directory = "{outdir}/{variety}.star.index",
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: 
        bam = "{outdir}/{sample}.{variety}.star.Aligned.sortedByCoord.out.bam",
    shell: r"""
    STAR \
    --readFilesIn {input.read_1} {input.read_2} \
    --readFilesCommand zcat \
    --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star. \
    --outSAMtype BAM SortedByCoordinate \
    --limitGenomeGenerateRAM 100000000000 \
    --limitBAMsortRAM 3000000000 \
    --genomeDir {input.star_index_directory}

    # Clean up uneccessary files 
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star.Log.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star.Log.final.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star.Log.progress.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star.SJ.out.tab
    """


rule star_rename:
    input: "{outdir}/{sample}.{variety}.star.Aligned.sortedByCoord.out.bam"
    output: "{outdir}/{sample}.{variety}.star.bam"
    shell: "mv {input} {output}"


rule bowtie2_align:
    # snakemake ../../../../../tmp/LE1.test.bowtie2.bam -s align.smk --cores 16
    input:
        bowtie2_index = multiext("{outdir}/{variety}.bowtie2.index.",            
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",),
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.bowtie2.bam"
    shell: r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.variety}.bowtie2.index -1 {input.read_1} -2 {input.read_2} | samtools view -bS > {output}
        """


rule bowtie2_align_large:
    input:
        bowtie2_index = multiext("{outdir}/{variety}.bowtie2.index.",            
            "1.bt2l",
            "2.bt2l",
            "3.bt2l",
            "4.bt2l",
            "rev.1.bt2l",
            "rev.2.bt2l",),
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.bowtie2.large.bam"
    shell: r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.variety}.bowtie2.index -1 {input.read_1} -2 {input.read_2} | samtools view -bS > {output}
        """


rule bwa_align:
    input:
        bwa_index = multiext("{outdir}/{variety}.bwa.index.",
        "amb",
        "ann",
        "bwt",
        "pac",
        "sa",),
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.bwa.bam"
    shell: r""" 
    bwa mem {wildcards.outdir}/{wildcards.variety}.bwa.index {input.read_1} {input.read_2} | samtools view -bS > {output}
    """


rule hisat2_align:
    input: 
        index = expand("{{outdir}}/{{variety}}.hisat2.index.{num}.ht2", num=["1", "2", "3", "4", "5", "6", "7", "8"]),
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.hisat2.bam"
    shell: r"""
    hisat2 -x {wildcards.outdir}/{wildcards.variety}.hisat2.index -p {threads} -1 {input.read_1} -2 {input.read_1} | \
    samtools view -bS > {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.hisat2.bam
    """


rule bamtools_statistics:
    input: "{outdir}/{sample}.{variety}.{tool}.bam"
    output: "{outdir}/{sample}.{variety}.{tool}.statistics.txt"
    shell: r"bamtools stats -in {input} > {output}"

rule report_stats_csv:
    input: "{outdir}/{sample}.{variety}.{tool}.statistics.txt"
    output: "{outdir}/{sample}.{variety}.{tool}.statistics.csv"
    script:
        "../scripts/alignment_stats.py"



#rule df_of_alignment_statistics:
#    input: expand("{{outdir}}/{sample}{repeats}.{variety}.{tool}.statistics.txt", sample=["LE", "SP", "RO", "IS"],
#     repeats=["1", "2", "3"], variety=["test", "fielder", "CS"], tool=["star", "bowtie2", "bwa", "hisat2"])
#    output: "{outdir}/alignment.csv"
#    shell:r"""
#    # For each LE1.test.star.statistics.txt file pull out the 
#    """

