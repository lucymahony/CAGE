import pandas as pd
import matplotlib.pyplot as plt

def extract_variables_from_file(file_path):
    file_name = file_path.rsplit('/', 1)[-1]
    sample = file_name[:2]
    repeat = file_name[2]
    variety = file_name.split('.')[1]
    tool = file_name.split('.')[2]
    print(f'the sample is {sample} ,the repeat is {repeat}, the variety is {variety}, the tool is {tool}')
    return sample, repeat, variety, tool 
     

def main(snakemake):
    list_input_file_paths = snakemake.input
    print(list_input_file_paths)
    dfs = []
    for file_path in list_input_file_paths:
        sample, repeat, variety, tool = extract_variables_from_file(file_path)
        df = pd.read_csv(file_path)
        df['Sample'] = sample
        df['Repeat'] = repeat
        df['Tool'] = tool
        dfs.append(df)
    concatenated_df = pd.concat(dfs, ignore_index=True)
    concatenated_df.to_csv(snakemake.output.csv, index=False)


if __name__ == "__main__":
    main(snakemake)

