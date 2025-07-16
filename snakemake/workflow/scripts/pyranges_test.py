import pandas as pd
import pyranges as pr

# I first need to double check how pyranges works, particularly .join() and .intersect() methods.
test_gene_annotation = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [10000, 50, 100],
 "End": [11000, 59, 111], "ID": ["a", "b", "c"], "Strand": ["+", "-", "+"]})
print(f"test_gene_annotation is {test_gene_annotation}")
test_cage = pr.from_dict({"Chromosome": ["chr1"] * 3, "Start": [10000, 58, 57], "End": [10001, 59, 59], "Strand": ["+", "-", "+"]})

join = test_gene_annotation.join(test_cage, strandedness="same")
print(join)
# length of unique ID's in intersect
print(f'The number is =  {len(set(join.df["ID"]))}')




# Read in Tag Cluster as PyRanges. 
#list_of_cage_bed_files = ["../../../intermediate_data/snakemake_intermediate_data/shared_RO_fielder_at_least_2_overlaps.bed", 
#"../../../intermediate_data/snakemake_intermediate_data/shared_LE_fielder_at_least_2_overlaps.bed",
#"../../../intermediate_data/snakemake_intermediate_data/shared_SP_fielder_at_least_2_overlaps.bed",
#"../../../intermediate_data/snakemake_intermediate_data/shared_IS_fielder_at_least_2_overlaps.bed",
#"../../../intermediate_data/snakemake_intermediate_data/shared_R_cadenza_at_least_2_overlaps.bed",
#"../../../intermediate_data/snakemake_intermediate_data/shared_L_cadenza_at_least_2_overlaps.bed"]
#
list_of_cage_bed_files = ["../../../intermediate_data/snakemake_intermediate_data/shared_RO_fielder_at_least_2_overlaps.bed"]


def read_cage_files_as_pyranges(list_of_cage_bed_files):
    all_cage_dfs = []
    for cage_bed_file in list_of_cage_bed_files:
        cage_columns = ["Chromosome", "Start", "End", "Name", "Scores", "Strand", "TC_Starts", "TC_Ends"]
        cage_df = pd.read_csv(cage_bed_file, sep="\t", header=None, names=cage_columns, index_col=False)
        cage_df["Start"] = cage_df["Start"].astype(int)
        cage_df["End"] = cage_df["End"].astype(int)
        all_cage_dfs.append(cage_df)
    concatenated_cage_df = pd.concat(all_cage_dfs, ignore_index=True)
    cage_ranges = pr.PyRanges(concatenated_cage_df)
    return cage_ranges

cage_ranges = read_cage_files_as_pyranges(list_of_cage_bed_files)

def read_in_annotation_as_pyranges(gff_file):
    gff_columns = ["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes"]
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, names=gff_columns, comment="#")
    gff_df["Start"] = gff_df["Start"].astype(int)
    gff_df["End"] = gff_df["End"].astype(int)
    print(gff_df.head())
    # Filter annotation for transcripts
    gff_transcripts = gff_df[gff_df["Feature"] == "mRNA"].copy()
    return gff_transcripts

gff_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/data/wheat_genome/CS_v2.1/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
gff_transcripts = read_in_annotation_as_pyranges(gff_file)


# Pyranges of TSS
def extract_start(row):
    # Extract the start position of transcripts based on strand
    return row["End"] if row["Strand"] == "-" else row["Start"]

gff_transcripts["Transcript_Start"] = gff_transcripts.apply(extract_start, axis=1)
transcript_starts = gff_transcripts[["Chromosome", "Transcript_Start", "Strand", "Attributes"]]
transcript_starts.rename(columns={"Transcript_Start": "Start"}, inplace=True)
transcript_starts.loc[:,"End"] = transcript_starts["Start"] + 1  # Make it a 1bp interval
tss_ranges = pr.PyRanges(transcript_starts) # Convert to PyRanges object

# Pyranges of 'transcriptional unit'
def extract_transcriptional_unit_start(row):
    return row["Start"] if row["Strand"] == "-" else row["Start"] - 500
def extract_transcriptional_unit_end(row):
    return row["End"] if row["Strand"] == "-" else row["End"] + 500

gff_transcripts["Transcriptional_Unit_Start"] = gff_transcripts.apply(extract_transcriptional_unit_start, axis=1)
gff_transcripts['Transcriptional_Unit_End'] = gff_transcripts.apply(extract_transcriptional_unit_end, axis=1)
transcriptional_unit = gff_transcripts[["Chromosome", "Transcriptional_Unit_Start", "Transcriptional_Unit_End", "Strand", "Attributes"]]
transcriptional_unit.rename(columns={"Transcriptional_Unit_Start": "Start", "Transcriptional_Unit_End": "End"}, inplace=True)

transcriptional_unit_pr = pr.PyRanges(transcriptional_unit)


# Expressed genes are the overlap between transcriptional unit and cage
expressed = cage_ranges.join(transcriptional_unit_pr, how="first", strandedness="same")
expressed_genes = expressed.df['Attributes'].tolist()
expressed_genes = [gene.split(';', 1)[0].split('=')[1] for gene in expressed_genes]
print(f'The number of unique expressed genes is {len(set(expressed_genes))} e.g. has a cage tag overlapping with 500 upstream to the end')



# Perform intersection to find overlaps between the CAGE ranges and the TSS.

overlaps = cage_ranges.join(tss_ranges, how="first", strandedness="same")
# Trun the Attributes column into a list and get in the  string in the " "
genes = overlaps.df['Attributes'].tolist()
genes = [gene.split(';', 1)[0].split('=')[1] for gene in genes]
unique_genes = len(set(genes))
print(f'Number of unique genes: {unique_genes} e.g. has a cage tag overlapping with the TSS')

# Number where the dominate tag is in the TSS


def read_cage_files_as_pyranges_dominant(list_of_cage_bed_files):
    all_cage_dfs = []
    for cage_bed_file in list_of_cage_bed_files:
        cage_columns = ["Chromosome", "Start", "End", "Name", "Scores", "Strand", "TC_Starts", "TC_Ends"]
        cage_df = pd.read_csv(cage_bed_file, sep="\t", header=None, names=cage_columns, index_col=False)
        cage_df["Start"] = cage_df["Start"].astype(int)
        cage_df["End"] = cage_df["End"].astype(int)
        all_cage_dfs.append(cage_df)
    concatenated_cage_df = pd.concat(all_cage_dfs, ignore_index=True)

    def filter_tc_starts_by_score(row):
        scores = list(map(float, row['Scores'].split(',')))
        tc_starts = row['TC_Starts'].split(',')
        max_score_index = scores.index(max(scores))
        return tc_starts[max_score_index]
    
    concatenated_cage_df["Dominant_TC_Start"] = concatenated_cage_df.apply(filter_tc_starts_by_score, axis=1)
    concatenated_cage_df.drop(columns=["Scores", "Start", "End", "TC_Starts", "TC_Ends"], inplace=True)
    print('Concatenated Cage DF with dominant tag cluster start')
    print(concatenated_cage_df.head())
    concatenated_cage_df.rename(columns={"Dominant_TC_Start": "Start"}, inplace=True)
    concatenated_cage_df["Start"] = concatenated_cage_df["Start"].astype(int)
    concatenated_cage_df.loc[:,"End"] = concatenated_cage_df["Start"] + 1  # Make it a 1bp interval
    dominant_tss = pr.PyRanges(transcript_starts) # Convert to PyRanges object
    return dominant_tss


cage_dominant_tss = read_cage_files_as_pyranges_dominant(list_of_cage_bed_files)
overlaps = cage_dominant_tss.join(tss_ranges, how="first", strandedness="same")
print('Dominant overlaps with TSS')
print(overlaps.df.head())
print(overlaps.df['Attributes'][:10])
# Trun the Attributes column into a list and get in the  string in the " "
genes = overlaps.df['Attributes'].tolist()
genes = [gene.split(';', 1)[0].split('=')[1] for gene in genes]
print(genes[:10])
unique_genes = len(set(genes))

print(f'Number of unique genes: {unique_genes} e.g. has a Dominant TSS tag overlapping with the TSS')
