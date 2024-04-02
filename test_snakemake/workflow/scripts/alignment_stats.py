# This file takes the output of bamtools statistics and turns it into a csv so that the results can be plotted 
import pandas as pd
import matplotlib.pyplot as plt 
import sys

metrics = ["Total reads", "Mapped reads", "Forward strand", "Reverse strand", "Failed QC", "Duplicates", "Paired-end reads"]
sys.stderr = open(snakemake.log[0], "w")

# Get the read-length
values = []
with open(snakemake.input["read_length"], 'r') as file:
    for line in file:
        for key in metrics:
            if line.startswith(key):
                value = line.split(":")[1].split()[0]
                values.append(value)

df = pd.DataFrame(values)
df.to_csv(snakemake.output["file_path"], index=False)



