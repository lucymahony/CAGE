
rule sort_bam:
    # Rule to sort BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.bam"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.bam")
    shell:
        "samtools sort  {input} -o {output}"

rule bam_to_sam:
    # Rule to convert sorted BAM files to SAM
    input:
    output:
        sorted_sam="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.sam"
    shell:
        "samtools view -h {input.sorted_bam} > {output.sorted_sam}"

rule generate_statistics:
    # Rule to generate summary statistics from BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.bam"
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
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.bam"
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
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.merge.ctss.bed")
    shell:
        "cat {input.plus_ctss} {input.minus_ctss} > {output}"

rule sort_merged_bed:
    # Rule to sort the merged ctss bed file
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.merge.ctss.bed"
    output: temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.bed")
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


rule standard_notation:
    # Rule to convert scientific notation to standard notation
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.bed"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.n.bed"
    shell:
        """
        # Use awk to process the file in a single pass:
        # 1. Match lines containing 'e' (scientific notation).
        # 2. Convert the 4th column (or whichever column contains scientific notation) to standard notation.
        # 3. Print modified lines with standard notation and unmodified lines directly.
        # 4. Add a header to the file
        
        awk '
        BEGIN {{
            OFS = "\\t"
        }}
        {{
            # Check if the fourth column contains scientific notation (contains "e").
            if ($4 ~ /e/) {{
                # Convert the scientific notation to standard format and update the fourth column.
                $4 = sprintf("%.0f", $4)
            }}
            # Print the updated line (modified or unmodified).
            print
        }}
        ' {input} > {output}
        """

rule add_headings:
    # Rule also adds a header to the file 'seqnames start strand count'. This is no longer required for the r script.
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.n.bed"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.h.bed"
    shell:
        """ 
        echo "seqnames    start    strand    count" > header.txt
        cat header.txt {input} > {output}
        rm header.txt
        """


rule update_statistics:
    # Rule to update the summary statistics with additional information
    # sorted_bed_file has the .n. to ensure that the file is in standard notation, which is required for downstream analysis with the R script.
    input:
        summary_statistics="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics.txt",
        plus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.ctss.bed",
        minus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.ctss.bed",
        merged_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.ctss.bed",
        sorted_bed_file="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.n.bed",
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
