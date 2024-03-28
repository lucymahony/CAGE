from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path

##### load config and sample sheets #####
samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index("sample", drop=False)

# for now I have commented out the actual samples 
# the test files aren't commented out and 
configfile: "config/config.yaml"


##### helper functions  #####
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_single_end(wildcards.sample, wildcards.unit):
        return units.loc[(wildcards.sample, wildcards.unit), "fq1"]
    else:
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        return [f"{u.fq1}", f"{u.fq2}"]


def get_untrim_counts(wildcards):
    """Get the counts for the untrimmed reads"""
    return f"{wildcards.tissue}_{read_number}.untrim.count"

def get_trim_counts(wildcards):
    """Get the counts for the trimmed reads"""
    return f"{wildcards.tissue}_{read_number}.trim.count"

