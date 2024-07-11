### Aligning the reads, first create the tools index files and then align them 
# Possibly the assemby and annotation names should be saved as tuples in the config file?
# For fielder the assemly name is  ____ for fielder the annotation is ________- note I need to add the liftover into the pipeline?
# Index files are called assembly_name.indexer_tool_name.index (maybe using the annotation name would be equivalent?)
# The test files are a small bacterial genome, just used in debugging test files. 


# At some point move these helper functions to the config file?
def get_fasta(genome_aligned_to):
    variety_assembly_dict = {'CS':'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna', # Note modified to have the same chromosome names as the annotation
                    'fielder': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/fielder_genome_data/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta', 
                    'cadenza': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/cadenza_genome_data/'} # Note currently no cadenza genome files
    return variety_assembly_dict[genome_aligned_to]

def get_gtf(genome_aligned_to):
    variety_annotation_dict = {'CS': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3', 
                    'fielder': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/fielder_genome_data/fielder.release.gtf',
                    'cadenza': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/cadenza_genome_data/'} # Note currently no cadenza genome files
    return variety_annotation_dict[genome_aligned_to]


########## Genome Indexing ##########

rule star_genome_index:
    # Generates a star index of given genome. 
    # NOTE that the names of the chromosomes on the assembly and the annotation have to be the same 
    # So for CS the assembly is modified using sed e.g. sed -e 's/>NC_057794.1 />Chr1A /' GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna
    # sed -e 's/>NC_057794.1 />Chr1A /
    
    input:
        fasta = lambda wildcards : get_fasta(wildcards.genome_aligned_to),  
        gtf = lambda wildcards : get_gtf(wildcards.genome_aligned_to),
    output: directory("{outdir}/{genome_aligned_to}.star.index")
    threads: 16
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
            --limitGenomeGenerateRAM 172020153610 \
            --genomeDir $tempdir \
        
        # Move outputs to proper location
        mv $tempdir {output} 
        """


rule bowtie2_genome_index:
    input: lambda wildcards : get_fasta(wildcards.genome_aligned_to)
    output: 
        multiext(
            "{outdir}/{genome_aligned_to}.bowtie2.index.",
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",
         )
    shell: 
        r"""
        bowtie2-build {input} {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index
        """


rule bowtie2_genome_index_large:
    input: lambda wildcards : get_fasta(wildcards.genome_aligned_to)
    output: 
        multiext(
            "{outdir}/{genome_aligned_to}.bowtie2.index.",
            "1.bt2l",
            "2.bt2l",
            "3.bt2l",
            "4.bt2l",
            "rev.1.bt2l",
            "rev.2.bt2l",
         )
    shell: 
        r"""
        bowtie2-build {input} {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index
        """


rule bwa_genome_index:
    input: lambda wildcards : get_fasta(wildcards.genome_aligned_to)
    output: 
        multiext(
            "{outdir}/{genome_aligned_to}.bwa.index.",
            "amb",
            "ann",
            "bwt",
            "pac",
            "sa"), 
    shell:
        r"""
        bwa index {input} -p {wildcards.outdir}/{wildcards.genome_aligned_to}.bwa.index
        """


rule hisat2_genome_index:
    input:
        fasta = lambda wildcards : get_fasta(wildcards.genome_aligned_to),  
        gtf = lambda wildcards : get_gtf(wildcards.genome_aligned_to),
    output: expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=[1, 2, 4, 5, 6, 7, 8])
    shell:
        r""" 
        hisat2_extract_splice_sites.py {input.gtf} > {wildcards.outdir}/fielder.ss
        hisat2_extract_exons.py {input.gtf} > {wildcards.outdir}/fielder.exon
        hisat2-build --ss {wildcards.outdir}/fielder.ss --exon {wildcards.outdir}/fielder.exon {input.fasta} {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index
        """


########## Aligning Paired Trimmed Reads  ##########

rule star_align_pe:
    # Aligning paired end reads with STAR
    input:
        star_index_directory = "{outdir}/{genome_aligned_to}.star.index",
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: 
        bam = "{outdir}/{sample}.{genome_aligned_to}.pe.star.Aligned.sortedByCoord.out.bam",
    shell: r"""
    STAR \
    --readFilesIn {input.read_1} {input.read_2} \
    --readFilesCommand zcat \
    --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star. \
    --outSAMtype BAM SortedByCoordinate \
    --limitGenomeGenerateRAM 100000000000 \
    --limitBAMsortRAM 3000000000 \
    --genomeDir {input.star_index_directory}

    # Clean up uneccessary files 
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.final.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.progress.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.SJ.out.tab
    """

rule star_align_se:
    # Aligning single-end reads with STAR
    input:
        star_index_directory = "{outdir}/{genome_aligned_to}.star.index",
        read = "{outdir}/{sample}.se.trim.fastq.gz",
    output: 
        bam = "{outdir}/{sample}.{genome_aligned_to}.se.star.Aligned.sortedByCoord.out.bam",
    shell: r"""
    STAR \
    --readFilesIn {input.read} \
    --readFilesCommand zcat \
    --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star. \
    --outSAMtype BAM SortedByCoordinate \
    --limitGenomeGenerateRAM 100000000000 \
    --limitBAMsortRAM 3000000000 \
    --genomeDir {input.star_index_directory}

    # Clean up unnecessary files 
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.final.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.Log.progress.out
    rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.star.SJ.out.tab
    """



rule star_rename:
    # This rule cleans up the output of star containing .Aligned.sortedByCoord.out in the output 
    # .{read_type} is either .pe or .se
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.star.Aligned.sortedByCoord.out.bam"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.star.bam"
    shell: "mv {input} {output}"


rule bowtie2_align_pe:
    # Aligning paired-end reads with bowtie 
    input:
        bowtie2_index = multiext("{outdir}/{genome_aligned_to}.bowtie2.index.",            
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",),
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: "{outdir}/{sample}.{genome_aligned_to}.pe.bowtie2.bam"
    shell: r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -1 {input.read_1} -2 {input.read_2} | samtools view -bS > {output}
        """

rule bowtie2_align_se:
    # Aligning single-end reads with bowtie 
    input:
        bowtie2_index = multiext("{outdir}/{genome_aligned_to}.bowtie2.index.",            
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",),
        read = "{outdir}/{sample}.se.trim.fastq.gz", 
    output: "{outdir}/{sample}.{genome_aligned_to}.se.bowtie2.bam"
    shell: r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -U {input.read} | samtools view -bS > {output}
        """

rule bowtie2_align_large:
    input:
        bowtie2_index = multiext("{outdir}/{genome_aligned_to}.bowtie2.index.",            
            "1.bt2l",
            "2.bt2l",
            "3.bt2l",
            "4.bt2l",
            "rev.1.bt2l",
            "rev.2.bt2l",),
        read_1 = "{outdir}/{sample}_1P.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.trim.fastq.gz",
    output: "{outdir}/{sample}.{genome_aligned_to}.bowtie2.large.bam"
    shell: r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -1 {input.read_1} -2 {input.read_2} | samtools view -bS > {output}
        """


rule bwa_align_pe:
    # Aligning paired-end reads with bwa
    input:
        bwa_index = multiext("{outdir}/{genome_aligned_to}.bwa.index.",
        "amb",
        "ann",
        "bwt",
        "pac",
        "sa",),
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: "{outdir}/{sample}.{genome_aligned_to}.pe.bwa.bam"
    shell: r""" 
    bwa mem {wildcards.outdir}/{wildcards.genome_aligned_to}.bwa.index {input.read_1} {input.read_2} | samtools view -bS > {output}
    """

rule bwa_align_se:
    # Aligning single-end reads with bwa
    input:
        bwa_index = multiext("{outdir}/{genome_aligned_to}.bwa.index.",
        "amb",
        "ann",
        "bwt",
        "pac",
        "sa",),
        read = "{outdir}/{sample}.se.trim.fastq.gz", 
    output: "{outdir}/{sample}.{genome_aligned_to}.se.bwa.bam"
    shell: r""" 
    bwa mem {wildcards.outdir}/{wildcards.genome_aligned_to}.bwa.index {input.read} | samtools view -bS > {output}
    """


rule hisat2_align_pe:
    # Align paired-end reads with hisat2 
    input: 
        index = expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=["1", "2", "3", "4", "5", "6", "7", "8"]),
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: "{outdir}/{sample}.{genome_aligned_to}.pe.hisat2.bam"
    shell: r"""
    hisat2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index -p {threads} -1 {input.read_1} -2 {input.read_1} | \
    samtools view -bS > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.hisat2.bam
    """


rule hisat2_align_se:
    # Align single-end reads with hisat2 
    input: 
        index = expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=["1", "2", "3", "4", "5", "6", "7", "8"]),
        read = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
    output: "{outdir}/{sample}.{genome_aligned_to}.pe.hisat2.bam"
    shell: r"""
    hisat2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index -p {threads} -U {input.read} | \
    samtools view -bS > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.hisat2.bam
    """


rule bamtools_statistics:
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.bam"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.statistics.txt"
    shell: r"bamtools stats -in {input} > {output}"


rule report_stats_csv:
    input: 
        statistics="{outdir}/{sample}.{genome_aligned_to}.{tool}.statistics.txt"
    output: 
        csv="{outdir}/{sample}.{genome_aligned_to}.{tool}.statistics.csv"
    script:
        "../scripts/alignment_stats.py"


rule df_of_alignment_statistics:
    # snakemake -j1 ../../../../../tmp/summary_alignment_stats_for_test.csv -s align.smk --use-conda --conda-frontend conda
    input: 
        #expand("{{outdir}}/{samples}{repeats}.{{genome_aligned_to}}.{tools}.statistics.csv", samples=["LE", "SP", "RO", "IS"], repeats=["1", "2", "3"], tools=["star", "bowtie2", "bwa", "hisat2"])
        expand("{{outdir}}/{samples}{repeats}.{{genome_aligned_to}}.{tools}.statistics.csv", samples=["LE"], repeats=["1"], tools=["star", "bowtie2"])
    output:
        csv="{outdir}/summary_alignment_stats_for_{genome_aligned_to}.csv"
    script:
        "../scripts/all_alignment_stats.py"


rule generate_alignment_plots:
    input:
        csv="{outdir}/summary_alignment_stats_for_{genome_aligned_to}.csv"
    output:
        "{outdir}/alignment_plot_{genome_aligned_to}.pdf"
    conda:
        "../envs/rscript.yaml"
    script:
        "../scripts/plot_alignment_results.R"
    
    
    
    
