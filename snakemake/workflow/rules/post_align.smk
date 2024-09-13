
rule sort_bam:
    # Rule to sort BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.bam"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_sorted.bam")
    shell:
        "samtools sort --threads 1 {input} -o {output}"

rule bam_to_sam:
    # Rule to convert sorted BAM files to SAM
    input:
        sorted_bam="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_sorted.bam"
    output:
        sorted_sam="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_sorted.sam"
    shell:
        "samtools view -h --threads 1 {input.sorted_bam} > {output.sorted_sam}"

rule generate_statistics:
    # Rule to generate summary statistics from BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_sorted.bam"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics.txt")
    shell:
        r"""
        echo "Total number of reads in the BAM file" > {output}
        samtools view -c {input} >> {output}

        echo "Counting only mapped (primary aligned) reads" >> {output}
        samtools view -c -F 260 {input} >> {output}
        """

rule bedtools_genomecov:
    # Rule to generate bed files for each strand using bedtools
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_sorted.bam"
    output:
        plus_bed=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.bed"),
        minus_bed=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.bed")
    shell:
        r"""
        bedtools genomecov -ibam {input} -5 -strand + -bg > {output.plus_bed}
        bedtools genomecov -ibam {input} -5 -strand - -bg > {output.minus_bed}
        """

rule ctss_conversion:
    # Rule to convert bed files to ctss format
    input:
        plus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.bed",
        minus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.bed"
    output:
        plus_ctss=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.ctss.bed"),
        minus_ctss=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.ctss.bed")
    shell:
        r"""
        awk '{{print $1 "\t" $2 "\t" "+" "\t" $4}}' {input.plus_bed} > {output.plus_ctss}
        awk '{{print $1 "\t" $2 "\t" "-" "\t" $4}}' {input.minus_bed} > {output.minus_ctss}
        """

rule merge_ctss:
    # Rule to merge ctss bed files
    input:
        plus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.ctss.bed",
        minus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.ctss.bed"
    output:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.ctss.bed"
    shell:
        "cat {input.plus_ctss} {input.minus_ctss} > {output}"

rule sort_merged_bed:
    # Rule to sort the merged ctss bed file
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.ctss.bed"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.bed"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"

rule update_statistics:
    # Rule to update the summary statistics with additional information
    input:
        summary_statistics="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics.txt",
        plus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.ctss.bed",
        minus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.ctss.bed",
        merged_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.ctss.bed",
        sorted_bed_file="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.bed",
        plus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.bed",
        minus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.bed",

    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics_updated.txt"
    shell:
        r"""
        cp {input.summary_statistics} {output} 

        echo "Statistics for {input.plus_bed}:" >> {output}
        wc -l {input.plus_bed} >> {output}

        echo "Statistics for {input.minus_bed}:" >> {output}
        wc -l {input.minus_bed} >> {output}

        echo "Statistics for {input.plus_ctss}:" >> {output}
        wc -l {input.plus_ctss} >> {output}

        echo "Statistics for {input.minus_ctss}:" >> {output}
        wc -l {input.minus_ctss} >> {output}

        echo "Statistics for {input.merged_ctss}:" >> {output}
        wc -l {input.merged_ctss} >> {output}

        echo "Statistics for {input.sorted_bed_file}:" >> {output}
        wc -l {input.sorted_bed_file} >> {output}

        echo "---------------------------------------------" >> {output}
        echo "" >> {output}
        """

rule conncatenate_statistics:
    input: getAllPostAlignmentStatistics
    output: "{outdir}/{genome_aligned_to}.{mapping_tool}.all_post_align_summary_statistics.txt"
    shell:
        "cat {input} > {output}"
