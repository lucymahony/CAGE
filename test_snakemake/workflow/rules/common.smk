from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path
from pandas.core.common import flatten
import os 

samples = pd.read_csv('../../config/samples.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)
units = pd.read_csv('../../config/units.tsv',  delim_whitespace=True, dtype=str, comment="#", index_col=False)

configfile: "../../config/config.yaml"

def get_fastqs(wildcards):
    """
    Get raw FASTQ files from unit sheet.
    """
    sample = wildcards.sample
    number_of_reads = units[units['sample'] == sample]['number_reads'].values[0] 
    print(f"The row in the units tsv is {units[units['sample'] == sample]}")


    if number_of_reads == '1': # Reads are single ended. 
        read_path = units[units['sample'] == sample]['fq1'].values[0]
        print(f'The read path is {read_path}')
        return [read_path]
    elif number_of_reads =='2': # The reads are paired ended 
        read_1, read_2 = units[units['sample'] == sample]['fq1'].values[0], units[units['sample'] == sample]['fq2'].values[0]
        print(f'The read paths are {read_1, read_2}')
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
    return untrimmed
