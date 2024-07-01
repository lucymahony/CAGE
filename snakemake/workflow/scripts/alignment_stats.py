# This file takes the output of bamtools statistics and turns it into a csv so that the results can be plotted 
import pandas as pd
import matplotlib.pyplot as plt 


def main(snakemake):
    list_of_metrics = ["Total reads", "Mapped reads", "Forward strand", "Reverse strand", "Failed QC", "Duplicates", "Paired-end reads"]
    input_file_path = snakemake.input.statistics
    dictionary = {}
    for metric in list_of_metrics:
        with open(input_file_path, 'r') as f:
            for line in f:
                if line.startswith(metric):

                    result = line.split(":")[1] # Take text after : (e.g. Mapped reads:      0    (0%) becomes       0    (0%))
                    result = result.lstrip() # (becomes 0    (0%))
                    result = int(result.split()[0]) # Take before white space (becomes 0)

                    dictionary[metric] = result

    df = pd.DataFrame(dictionary, index=[0])
    output_file_path = snakemake.output.csv
    df.to_csv(output_file_path, index=False)


if __name__ == "__main__":
    main(snakemake)

