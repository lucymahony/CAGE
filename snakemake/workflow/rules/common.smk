from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path
from pandas.core.common import flatten
import os 

samples = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/samples.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)
units = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/units.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)

configfile: "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/snakemake/config/config.yaml"

def get_fastqs(wildcards):
    """
    Get raw FASTQ files from unit sheet.
    """
    sample = wildcards.sample.split('_R')[0] # Removes the bit in the case of LE1_R1 e.c.t 
    matching_units = units[units['sample'] == sample]
    if matching_units.empty:
        raise ValueError(f"No entry found for sample {sample} in units.tsv")
        
    number_of_reads = units[units['sample'] == sample]['number_reads'].values[0] 
    if number_of_reads == '1': # Reads are single ended. 
        read_path = units[units['sample'] == sample]['fq1'].values[0]
        return [read_path]
    elif number_of_reads =='2': # The reads are paired ended 
        read_1, read_2 = units[units['sample'] == sample]['fq1'].values[0], units[units['sample'] == sample]['fq2'].values[0]
        return [read_1, read_2]
    else: print('Error getting the number of reads from the units.tsv file')


def getAllTrimmedCountFiles(wildcards):
    """
    Returns list of trimmed count files 
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()
    path = wildcards.outdir

    x = [[sample + '_R1.pe', sample + '_R2.pe'] if number == "2" else [sample + '.se'] for sample, number in zip(samples, reads)] 
    flat  = [item for sublist in x for item in sublist]
    trimmed = [os.path.join(path, name + '.trim.count') for name in flat] 
    print(f'all trimmed files are {trimmed}')
    return trimmed


def getAllUntrimmedCountFiles(wildcards):
    """
    Returns list of untrimmed count files 
    Snakemake input functions can only return single or lists of not dictionaries and thereore there is a 
    getAllTrimmedCountFiles and a getAllUntrimmedCountFiles function
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()
    path = wildcards.outdir
    x = [[sample + '_R1.pe', sample + '_R2.pe'] if number == "2" else [sample + '.se'] for sample, number in zip(samples, reads)] 
    flat  = [item for sublist in x for item in sublist]
    untrimmed = [os.path.join(path, name + '.untrim.count') for name in flat]
    print(f'all untrimmed files are {untrimmed}')
    return untrimmed


def getAllAllignmentStatisticCSVs(wildcards):
    """
    Returns a list of the .CS.se.star.statistics.csv from the samples in samples.tsv for a given genome aligned to and aligner tool, which are set in the output file name as wildcards
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()

    genome_aligned_to = wildcards.genome_aligned_to
    alignment_tool = wildcards.tool
    path = wildcards.outdir

    alignment_statistics = []

    for sample, number in zip(samples, reads):
        if number == '1':
            alignment_statistics.append(f"{sample}.{genome_aligned_to}.se.{alignment_tool}.statistics.csv")
        else:
            endings = ['R1', 'R2']
            alignment_statistics.extend([
                f"{sample}_{end}.{genome_aligned_to}.pe.{alignment_tool}.statistics.csv" for end in endings
            ])
    alignment_statistics = [os.path.join(path, name) for name in alignment_statistics]
    print(f'the alignment statistics are {alignment_statistics}')
    return alignment_statistics


def getAllPostAlignmentStatistics(wildcards):
    """
    Returns a list of the .post_align_summary_statistics.txt from the samples in samples.tsv for a given genome aligned to and aligner tool, which are set in the output file name as wildcards
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()
    print(f'the reads are {reads}')

    genome_aligned_to = wildcards.genome_aligned_to
    mapping_tool = wildcards.mapping_tool
    path = wildcards.outdir

    post_alignment_statistics = []
    for sample, read in zip(samples, reads):
        if read == '1':
            read_type = 'se'
        else:
            read_type = 'pe'
        name = f"{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.post_align_summary_statistics_updated.txt"
        post_alignment_statistics.append(os.path.join(path, name))
    print(f'the post alignment statistics are {post_alignment_statistics}')
    return post_alignment_statistics


def getAllSortedCTSSNBedfiles(wildcards):
    """
    Returns a list of the sorted CTSS bed files from the samples in samples.tsv for a given genome aligned to and aligner tool, which are set in the output file name as wildcards
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()

    genome_aligned_to = wildcards.genome_aligned_to
    mapping_tool = wildcards.mapping_tool
    path = wildcards.outdir

    sorted_ctss_bedfiles = []
    for sample, read in zip(samples, reads):
        if read == '1':
            read_type = 'se'
        else:
            read_type = 'pe'
        name = f"{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}.sorted.ctss.n.bed"
        sorted_ctss_bedfiles.append(os.path.join(path, name))
    return sorted_ctss_bedfiles


def getAllBigWigFiles(wildcards):
    """
    Returns a list of bigwig files for CAGE fight r analysis 
        
    """
    samples = units['sample'].tolist()
    reads = units['number_reads'].tolist()

    genome_aligned_to = wildcards.genome_aligned_to
    mapping_tool = wildcards.mapping_tool
    path = wildcards.outdir
    read_type = wildcards.read_type

    bw_plus_files = []
    bw_minus_files = []
    for sample, read in zip(samples, reads):
        name = f"{sample}.{genome_aligned_to}.{read_type}.{mapping_tool}"
        bw_plus_files.append(os.path.join(path, name) + '_plus.unique.bw')
        bw_minus_files.append(os.path.join(path, name) + '_minus.unique.bw')
    return bw_plus_files, bw_minus_files
