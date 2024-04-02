### Aligning the reads, first create the tools index files and then align them 
# Possibly the assemby and annotation names should be saved as tuples in the config file?
# For fielder the assemly name is  ____ for fielder the annotation is ________- note I need to add the liftover into the pipeline?
# Index files are called assembly_name.indexer_tool_name.index (maybe using the annotation name would be equivalent?)
# For chinese spring the assembly and annotation are called 
# test_CAGE_data/genomes/iwgsc_refseqv2.1_annotation_200916_HC.gff3 and  iwgsc_refseqv2.1_assembly.fa 


# At some point move these helper functions to the config file?
def get_fasta(variety):
    variety_assembly_dict = {'CS': '../../../test_CAGE_data/genomes/iwgsc_refseqv2.1_assembly.fa', 
                    'fielder': '../../../test_CAGE_data/genomes/fielder.fa'}
    return variety_assembly_dict[variety]

def get_gtf(variety):
    variety_annotation_dict = {'CS': '../../../test_CAGE_data/genomes/iwgsc_refseqv2.1_annotation_200916_HC.gff3', 
                    'fielder': '../../../test_CAGE_data/genomes/fielder.gtf'}
    return variety_annotation_dict[variety]


########## Genome Indexing ##########

rule star_genome_index:
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
            --genomeDir $tempdir
        
        # Move outputs to proper location
        mv $tempdir {output} 
        """


rule bowtie2_genome_index:
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
            "{outdir}/{variety}.bwa.index",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa"), 
    shell:
        r"""
        bwa index {input} -p {wildcards.outdir}/{wildcards.variety}.bwa.index
        """


rule hisat2_genome_index:
    input:
        fasta = lambda wildcards : get_fasta(wildcards.variety),  
        gtf = lambda wildcards : get_gtf(wildcards.variety),
    output: "{outdir}/{variety}.hisat2.index"
    shell:
        r""" 
        hisat2_extract_splice_sites.py {input.gtf} > {wildcards.outdir}/fielder.ss
        hisat2_extract_exons.py {input.gtf} > {wildcards.outdir}/fielder.exon
        hisat2-build --ss {wildcards.outdir}/fielder.ss --exon {wildcards.outdir}/fielder.exon {input.fasta} {output}
        """


########## Aligning Paired Trimmed Reads  ##########
#move the indexes into the tmp directory 

rule star_align:
    input:
        star_index_directory = "{outdir}/{variety}.star.index",
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.star.bam"

    shell: r"""
    STAR 
    --genomeDir {input.star_index_directory} \
    --readFilesIn {input.read_1} {input.read_2} \
    --readFilesCommand zcat \
    --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.variety}.star \
    --outSAMtype BAM SortedByCoordinate \
    --limitGenomeGenerateRAM 100000000000 \
    --limitBAMsortRAM 3000000000
    """


rule bowtie_align:
    input:
        bowtie_index = "{outdir}/bowtie_index_{genome}.idx",
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.bowtie.bam"
    shell: r""" 
        bowtie2 -p {threads} -x {input.bowtie_index} -1  {input.read_1} -2 {input.read_2} | samtools view -bS > {output}
        """


rule bwa_align:
    input:
        bwa_index = "{outdir}/bwa_index_{genome}.idx",
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.bwa.bam"
    shell: r""" 
    bwa mem {input.bwa_index} {input.read_1} {input.read_2} | samtools view -bS > {output}
    """


rule hisat_align:
    output: "{outdir}/hisat_{sample}_aligned.bam"
    input: 
        index = "{outdir}/hisat_index_{genome}.idx",
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{variety}.hisat.bam"
    shell: r"""
    hisat2 -x {input.index} -p {threads} -1 {input.read_1} -2 {input.read_1} | \
    samtools view -bS > {outdir}/{sample}.hisat.bam
    """

