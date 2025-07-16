
rule filter_uniquely_mapped_reads:
    input:
        bam="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.bam"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.bam")
    shell:
        """
        samtools view -h -q 255 {input.bam} > {output}
        """

rule sort_bam:
    # Rule to sort BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.bam"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.bam")
    shell:
        "samtools sort  {input} -o {output}"


rule generate_statistics:
    # Rule to generate summary statistics from BAM files
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.bam"
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
    input:"{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.bam"
    output:
        plus_bed=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.bed"),
        minus_bed=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.bed")
    
    shell:
        r"""
        bedtools genomecov -ibam {input} -5 -strand + -bg > {output.plus_bed}
        bedtools genomecov -ibam {input} -5 -strand - -bg > {output.minus_bed}
        """

rule create_genome_file:
    # Need a file of the genome sizes for bedGraphToBigWig - Taken from the forged genome 
    output:temp("{outdir}/Taestivum.ChineseSpring.chrom.sizes")
    shell:
        r"""
        # Create a file with the chromosome sizes
        echo -e "Chr1A\t598660471" > {output}
        echo -e "Chr1B\t700547350" >> {output}
        echo -e "Chr1D\t498638509" >> {output}
        echo -e "Chr2A\t787782082" >> {output}
        echo -e "Chr2B\t812755788" >> {output}
        echo -e "Chr2D\t725972874" >> {output}
        echo -e "Chr3A\t750843639" >> {output}
        echo -e "Chr3B\t830837813" >> {output}
        echo -e "Chr3D\t718151956" >> {output}
        echo -e "Chr4A\t744588211" >> {output}
        echo -e "Chr4B\t781547799" >> {output}
        echo -e "Chr4D\t650055668" >> {output}
        echo -e "Chr5A\t773760400" >> {output}
        echo -e "Chr5B\t805254273" >> {output}
        echo -e "Chr5D\t677843184" >> {output}
        echo -e "Chr6A\t731188232" >> {output}
        echo -e "Chr6B\t731188232" >> {output}
        echo -e "Chr6D\t495380293" >> {output}
        echo -e "Chr7A\t744491536" >> {output}
        echo -e "Chr7B\t764072961" >> {output}
        echo -e "Chr7D\t642921167" >> {output}
        """

rule bam_to_bigwig:
    # If using CAGEfightR downstream, having bigwig files is useful
    input:
        plus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.bed",
        minus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.bed",
        genome_file="{outdir}/Taestivum.ChineseSpring.chrom.sizes",
    output:
        plus_bigwig="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.bw",
        minus_bigwig="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.bw",
    shell:
        r"""
        # Filter BED files to remove unwanted chromosomes - Remove those not starting with Chr

        awk 'NR==FNR {{a[$1]; next}} $1 in a' {input.genome_file} {input.plus_bed} > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.plus.unique.bed.filtered
        awk 'NR==FNR {{a[$1]; next}} $1 in a'  {input.genome_file} {input.plus_bed} > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.minus.unique.bed.filtered

 
        

        # Generate coverage and convert to BigWig
        # Awk step is required so that width is consistently 1bp, e.g. Chr1A 5 7 1 -> Chr1A 5 6 1 Chr1A 6 7 1, which is required for CAGEfightR

        bedtools genomecov -i  {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.plus.unique.bed.filtered -bg -g {input.genome_file} > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_plus.unique.bed.bg
        awk '{{for(i=$2; i<$3; i++) print $1, i, i+1, $4}}' {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_plus.unique.bed.bg > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_plus.unique.bed.bg.tmp
        bedGraphToBigWig  {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_plus.unique.bed.bg.tmp {input.genome_file} {output.plus_bigwig}

        bedtools genomecov -i {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}.minus.unique.bed.filtered -bg -g {input.genome_file} > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_minus.unique.bed.bg
        awk '{{for(i=$2; i<$3; i++) print $1, i, i+1, $4}}' {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_minus.unique.bed.bg > {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_minus.unique.bed.bg.tmp
        bedGraphToBigWig  {wildcards.outdir}/{wildcards.sample}.{wildcards.genome_aligned_to}.{wildcards.read_type}.{wildcards.mapping_tool}_minus.unique.bed.bg.tmp {input.genome_file} {output.minus_bigwig}
        """


rule ctss_conversion:
    # Rule to convert bed files to ctss format
    input:
        plus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.bed",
        minus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.bed"
    output:
        plus_ctss=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.ctss.bed"),
        minus_ctss=temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.ctss.bed")
    shell:
        r"""
        awk '{{print $1 "\t" $2 "\t" "+" "\t" $4}}' {input.plus_bed} > {output.plus_ctss}
        awk '{{print $1 "\t" $2 "\t" "-" "\t" $4}}' {input.minus_bed} > {output.minus_ctss}
        """

rule merge_ctss:
    # Rule to merge ctss bed files
    input:
        plus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.ctss.bed",
        minus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.ctss.bed"
    output:temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.merge.ctss.bed")
    shell:
        "cat {input.plus_ctss} {input.minus_ctss} > {output}"

rule sort_merged_bed:
    # Rule to sort the merged ctss bed file
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.merge.ctss.bed"
    output: temp("{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.bed")
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


rule standard_notation:
    # Rule to convert scientific notation to standard notation
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.bed"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.n.bed"
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
    input: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.n.bed"
    output: "{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.h.bed"
    shell:
        """ 
        echo "seqnames    start    strand    count" > header.txt
        cat header.txt {input} > {output}
        """



rule update_statistics:
    # Rule to update the summary statistics with additional information
    # sorted_bed_file has the .n. to ensure that the file is in standard notation, which is required for downstream analysis with the R script.
    input:
        summary_statistics="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics.txt",
        plus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.ctss.bed",
        minus_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.ctss.bed",
        merged_ctss="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.merge.ctss.bed",
        sorted_bed_file="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.unique.sorted.ctss.n.bed",
        plus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_plus.unique.bed",
        minus_bed="{outdir}/{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}_minus.unique.bed",

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

