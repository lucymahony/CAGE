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

    x = [[sample + '_R1.pe', sample + '_R2.pe'] if number == "2" else [sample + '.se'] for sample, number in zip(samples, reads)] # E.g. IS1_R1
    flat  = [item for sublist in x for item in sublist]
    trimmed = [os.path.join(path, name + '.trim.count') for name in flat]  # E.g. ../../../../../tmp/IS1_R1.trim.count
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
    x = [[sample + '_R1.pe', sample + '_R2.pe'] if number == "2" else [sample + '.se'] for sample, number in zip(samples, reads)] # E.g. IS1_R1
    flat  = [item for sublist in x for item in sublist]
    untrimmed = [os.path.join(path, name + '.untrim.count') for name in flat]
    print(f'all untrimmed files are {untrimmed}')
    return untrimmed
