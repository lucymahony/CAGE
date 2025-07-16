# This script generates genome indexes, maps reads, and does post mapping processing for a variety of mapping tools.

########## Helper functions ########## 

def get_fasta(genome_aligned_to):
    variety_assembly_dict = {'CS':'/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna', # Note modified to have the same chromosome names as the annotation
                    'fielder': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/fielder_genome_data/201216_Fielder_pseudomolecules_V1+unanchored_contigs.fasta', 
                    'cadenza': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/cadenza_genome_data/24092024_cadenza_v2.fa'} 
    return variety_assembly_dict[genome_aligned_to]


def get_gtf(genome_aligned_to):
    variety_annotation_dict = {'CS': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3', 
                    'fielder': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/fielder_genome_data/fielder.release.gtf',
                    'cadenza': '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/cadenza_genome_data/TraesCAD_EIv1.0.release.gff3'} 
    return variety_annotation_dict[genome_aligned_to]



########## Genome Indexing ##########

rule star_genome_index:
    """
    Generates a star index of given genome. 
    The names of the chromosomes on the assembly and the annotation have to be the same 
    So for CS the assembly is modified using sed e.g. sed -e 's/>NC_057794.1 />Chr1A /' GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna
    sed -e 's/>NC_057794.1 />Chr1A /' \
        -e 's/>NC_057795.1 />Chr1B /' \
        -e 's/>NC_057796.1 />Chr1D /' \
        -e 's/>NC_057797.1 />Chr2A /' \
        -e 's/>NC_057798.1 />Chr2B /' \
        -e 's/>NC_057799.1 />Chr2D /' \
        -e 's/>NC_057800.1 />Chr3A /' \
        -e 's/>NC_057801.1 />Chr3B /' \
        -e 's/>NC_057802.1 />Chr3D /' \
        -e 's/>NC_057803.1 />Chr4A /' \
        -e 's/>NC_057804.1 />Chr4B /' \
        -e 's/>NC_057805.1 />Chr4D /' \
        -e 's/>NC_057806.1 />Chr5A /' \
        -e 's/>NC_057807.1 />Chr5B /' \
        -e 's/>NC_057808.1 />Chr5D /' \
        -e 's/>NC_057809.1 />Chr6A /' \
        -e 's/>NC_057810.1 />Chr6B /' \
        -e 's/>NC_057811.1 />Chr6D /' \
        -e 's/>NC_057812.1 />Chr7A /' \
        -e 's/>NC_057813.1 />Chr7B /' \
        -e 's/>NC_057814.1 />Chr7D /' \
        GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic_modified.fna 
        which is in this script bash CAGE/input_data/chinese_spring_genome_data/modify_assembly.sh
    """
    input:
        fasta = lambda wildcards : get_fasta(wildcards.genome_aligned_to),  
        gtf = lambda wildcards : get_gtf(wildcards.genome_aligned_to),
    output: directory("{outdir}/{genome_aligned_to}.star.index")
    threads: 16
    shell:
        r"""
        # Create temporary directory
        tempdir=$(mktemp -d)
        trap "rm -rf $tempdir" EXIT

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100 \
            --limitGenomeGenerateRAM 222020153610 \
            --genomeDir $tempdir \
        
        # Move outputs to proper location
        mv $tempdir {output} 
        """


rule bowtie2_genome_index:
    # Generates a bowtie2 index of given genome.
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
    threads: 16
    shell: 
        r"""
        bowtie2-build {input} {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index
        """


rule bowtie2_genome_index_large:
    # Generates a bowtie2 index of given genome, when the genome is large e.g. bt2l file ending rather than bt2
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
    # Generates a bwa index of a given genome
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
    # Generates a hisat2 index of a given genome
    input:
        fasta = lambda wildcards : get_fasta(wildcards.genome_aligned_to),  
        gtf = lambda wildcards : get_gtf(wildcards.genome_aligned_to),
    output: expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=[1, 2, 3, 4, 5, 6, 7, 8])
    shell:
        r""" 
        hisat2_extract_splice_sites.py {input.gtf} > {wildcards.outdir}/fielder.ss
        hisat2_extract_exons.py {input.gtf} > {wildcards.outdir}/fielder.exon
        hisat2-build --ss {wildcards.outdir}/fielder.ss --exon {wildcards.outdir}/fielder.exon {input.fasta} {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index
        """


########## Aligning Trimmed Reads ########## 

rule star_align_pe:
    # Aligning paired end reads with STAR
    input:
        star_index_directory = "{outdir}/{genome_aligned_to}.star.index",
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: 
        bam = "{outdir}/{sample}.{genome_aligned_to}.pe.star.Aligned.sortedByCoord.out.bam", # Check that this is correct _R1
    threads: 16
    shell: 
        r"""
        STAR \
        --runThreadN {threads} \
        --readFilesIn {input.read_1} {input.read_2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star. \
        --outSAMtype BAM SortedByCoordinate \
        --limitGenomeGenerateRAM 100000000000 \
        --limitBAMsortRAM 3000000000 \
        --genomeDir {input.star_index_directory}

        # Clean up uneccessary files 
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star.Log.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star.Log.final.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star.Log.progress.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star.SJ.out.tab
        rm -rf {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.star._STARtmp
        """


rule star_align_se:
    # Aligning single-end reads with STAR
    input:
        star_index_directory = "{outdir}/{genome_aligned_to}.star.index",
        read = "{outdir}/{sample}.se.trim.fastq.gz",
    output: 
        bam = "{outdir}/{sample}.{genome_aligned_to}.se.star.Aligned.sortedByCoord.out.bam",
    threads: 16
    shell: 
        r"""
        STAR \
        --runThreadN {threads} \
        --readFilesIn {input.read} \
        --readFilesCommand zcat \
        --outFileNamePrefix {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star. \
        --outSAMtype BAM SortedByCoordinate \
        --limitGenomeGenerateRAM 100000000000 \
        --limitBAMsortRAM 3000000000 \
        --genomeDir {input.star_index_directory}

        # Clean up unnecessary files 
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star.Log.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star.Log.final.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star.Log.progress.out
        rm {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star.SJ.out.tab
        rm -rf {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.star._STARtmp
        """


rule star_rename:
    # Cleans up the output name of the mapped reads from STAR 
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
    threads: 16
    shell: 
        r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -1 {input.read_1} -2 {input.read_2} --threads {threads}| samtools view -bS > {output}
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
    threads: 16
    shell: 
        r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -U {input.read} --threads {threads}| samtools view -bS > {output}
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
    threads: 16
    shell: 
        r""" 
        bowtie2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.bowtie2.index -1 {input.read_1} -2 {input.read_2} --threads {threads}| samtools view -bS > {output}
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
    threads: 16
    shell: 
        r""" 
        bwa mem {wildcards.outdir}/{wildcards.genome_aligned_to}.bwa.index {input.read_1} {input.read_2} -t {threads} | samtools view -bS > {output}
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
    threads: 16
    shell: 
        r""" 
        bwa mem {wildcards.outdir}/{wildcards.genome_aligned_to}.bwa.index {input.read} -t {threads} | samtools view -bS > {output}
        """


rule hisat2_align_pe:
    # Align paired-end reads with hisat2 
    input: 
        index = expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=["1", "2", "3", "4", "5", "6", "7", "8"]),
        read_1 = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
        read_2 = "{outdir}/{sample}_2P.pe.trim.fastq.gz",
    output: "{outdir}/{sample}.{genome_aligned_to}.pe.hisat2.bam"
    threads: 16
    shell: 
        r"""
        hisat2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index -p {threads} -1 {input.read_1} -2 {input.read_1} | \
        samtools view -bS > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.pe.hisat2.bam
        """


rule hisat2_align_se:
    # Align single-end reads with hisat2 
    input: 
        index = expand("{{outdir}}/{{genome_aligned_to}}.hisat2.index.{num}.ht2", num=["1", "2", "3", "4", "5", "6", "7", "8"]),
        read = "{outdir}/{sample}_1P.pe.trim.fastq.gz", 
    output: "{outdir}/{sample}.{genome_aligned_to}.se.hisat2.bam"
    threads: 16
    shell: 
        r"""
        hisat2 -x {wildcards.outdir}/{wildcards.genome_aligned_to}.hisat2.index -p {threads} -U {input.read} | \
        samtools view -bS > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.se.hisat2.bam
        """


rule bamtools_statistics:
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.bam"
    output: temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.statistics.txt") 
    shell: 
        r"bamtools stats -in {input} > {output}"


rule report_stats_csv:
    # Turns the text file output of bamtools statistics and turns it into a csv so that it can be concatenated in the rule df_of_alignment_statistics
    input: 
        statistics="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.statistics.txt"
    output: 
        csv=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{tool}.statistics.csv")
    script:
        "../scripts/alignment_stats.py"


rule df_of_alignment_statistics:
    # Generates a dataframe of the alignment statistics from bamtools of the files listed in units.csv
    input: 
        getAllAllignmentStatisticCSVs
    output:
        csv="{outdir}/summary_alignment_stats.{genome_aligned_to}.{tool}.csv"
    script:
        "../scripts/all_alignment_stats.py"


rule generate_alignment_plots:
    # Note, currently the r packages are installed in snakemake conda env. But can instead do conda: "../envs/rscript.yaml"
    input:
        csv="{outdir}/summary_alignment_stats_for_{genome_aligned_to}.csv"
    output:
        "{outdir}/alignment_plot_{genome_aligned_to}.pdf"
    script:
        "../scripts/plot_alignment_results.R"
        

