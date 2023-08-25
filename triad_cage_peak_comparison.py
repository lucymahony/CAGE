# This script takes a bed file  of CAGE Tag cluster peaks, a table of wheat triads, and a gff file with the wheat
# genome annotation
# and returns a dataframe that for each triad has a TRUE of FALSE if there is a CAGE peak within a given
# distance up and down stream from the gene start on the same strand

# Import libraries
import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np

# Import datasets
triad_table = pd.read_csv('/Users/mahony/Documents/1.wheat_triads_project/data/Triad_data/hc_triad_table copy',
                          sep='\t')
triad_table.columns = ['A', 'B', 'D']


def get_genome_annotation(chromosome):
    # genome_annotation = pd.read_csv(f'/Users/mahony/Documents/my_cage_tool_files/CS_annotation_files/'
    #                                 f'Triticum_aestivum.IWGSC.56.chromosome.{chromosome}.gff3', skiprows=6, sep='\t',
    #                                 comment='#')
    genome_annotation = pd.read_csv(f'/Users/mahony/Downloads/fielder.release.gff', sep='\t',
                                    comment='#')
    genome_annotation.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'ID']
    genome_annotation = genome_annotation[genome_annotation['type'] == 'gene']
    genome_annotation = genome_annotation[genome_annotation['chr'] == chromosome] # I've done this so that fielder and CS are treated the same
    # Strip ID=gene: from the ID column and split at the first ; to remove the rest of the annotation
    genome_annotation['ID'] = genome_annotation['ID'].str.replace('ID=', '') # ID=gene: for CS
    genome_annotation['ID'] = genome_annotation['ID'].str.split(';').apply(lambda x: x[0])
    return genome_annotation


def make_gene_dict(triad_table):
    """
    Function to make a dictionary of genomic information for each wheat gene in a triad. The dictionary is keyed by
    the gene name and the values are a list of lists of the genomic information for each gene
    chromosome, strand, gene_start
    """
    # Turn the triad table into a list of genes
    triad_genes = triad_table.values.flatten().tolist()
    # Empty dictionary
    gene_dict = {}
    # For efficiency read the genome annotation for each chromosome once
    # Create a set of unique chromosomes
    chromosomes = {gene[7:9] for gene in triad_genes}
    for chromosome in tqdm(chromosomes, desc=f'Getting genome annotation'):
        genome_annotation = get_genome_annotation(chromosome)
        # Get all triad genes on this chromosome
        triad_genes_on_chromosome = [gene for gene in triad_genes if gene[7:9] == chromosome]
        for gene in triad_genes_on_chromosome:
            gene_info = genome_annotation[genome_annotation['ID'] == gene]
            if gene_info.empty:
                continue  # Skip if gene not found in genome_annotation
            else:
                strand = gene_info['strand'].values[0]
                if strand == '+':
                    gene_start = gene_info['start'].values[0]
                else:
                    gene_start = gene_info['end'].values[0]
                gene_dict[gene] = [chromosome, strand, gene_start]
    # Return the completed dictionary
    return gene_dict


def make_gene_dict_from_list(gene_list):
    """
    Function to make a dictionary of genomic information for each wheat gene in a triad. The dictionary is keyed by
    the gene name and the values are a list of lists of the genomic information for each gene
    chromosome, strand, gene_start
    """
    # Turn the triad table into a list of genes

    # Empty dictionary
    gene_dict = {}
    # For efficiency read the genome annotation for each chromosome once
    # Create a set of unique chromosomes
    chromosomes = {gene[7:9] for gene in gene_list}
    for chromosome in tqdm(chromosomes, desc=f'Getting genome annotation'):
        genome_annotation = get_genome_annotation(chromosome)
        # Get all triad genes on this chromosome
        triad_genes_on_chromosome = [gene for gene in gene_list if gene[7:9] == chromosome]
        for gene in triad_genes_on_chromosome:
            gene_info = genome_annotation[genome_annotation['ID'] == gene]
            if gene_info.empty:
                continue  # Skip if gene not found in genome_annotation
            else:
                strand = gene_info['strand'].values[0]
                if strand == '+':
                    gene_start = gene_info['start'].values[0]
                else:
                    gene_start = gene_info['end'].values[0]
                gene_dict[gene] = [chromosome, strand, gene_start]
    # Return the completed dictionary
    return gene_dict


def check_for_cage_peak(chromosome, strand, region_start, region_stop, cage_peaks):
    """
    Function that checks if given region has a CAGE peak. NOTE AT THE MOMENT IM UNSURE ABOUT THE EDGE OF THE REGION
    IF IT SHOULD BE INCLUDED WITH <= OR NOT. CURRENTLY IT IS INCLUDED
    :param chromosome: chromosome of the region
    :param strand: strand of the region
    :param region_start:
    :param region_stop:
    :param cage_peaks: dataframe of cage peaks
    :return: TRUE or FALSE
    """
    # Filter the cage peaks to the chromosome and strand
    cage_peaks = cage_peaks[cage_peaks['seqnames'] == chromosome]
    cage_peaks = cage_peaks[cage_peaks['strand'] == strand]

    # Check overlaps of the cage start.ie the cage start is bigger than the region start but smaller than the region end
    subset = cage_peaks[(cage_peaks['start'] >= region_start) & (cage_peaks['start'] <= region_stop)]
    if not subset.empty:
        return True
    else:
        # Check overlaps at the cage end.ie the cage end is smaller than the region end but bigger than the region start
        subset = cage_peaks[(cage_peaks['end'] <= region_stop) & (cage_peaks['end'] >= region_start)]
        if not subset.empty:
            return True
        else:
            return False


def check_for_dominant_cage_peak(chromosome, strand, region_start, region_stop, cage_peaks):
    """
    Function same as check_for_cage_peak but checks for dominant peak rather than range
    """
    # Filter the cage peaks to the chromosome and strand
    cage_peaks = cage_peaks[cage_peaks['seqnames'] == chromosome]
    cage_peaks = cage_peaks[cage_peaks['strand'] == strand]

    # Check overlaps of the cage start.ie the cage start is bigger than the region start but smaller than the region end
    subset = cage_peaks[(cage_peaks['dominant_ctss'] >= region_start) & (cage_peaks['dominant_ctss'] <= region_stop)]
    if not subset.empty:
        return True
    else:
        # Check overlaps at the cage end.ie the cage end is smaller than the region end but bigger than the region start
        subset = cage_peaks[(cage_peaks['end'] <= region_stop) & (cage_peaks['end'] >= region_start)]
        if not subset.empty:
            return True
        else:
            return False


# Score transcript as True or False if there is a CAGE peak within a given distance of the gene.
# Get the chromosome, strand, gene_start from the dictionary
# turn Gene start into a region start and stop by adding and subtracting the distance
# Put through the check_for_cage_peak function

def score_transcript(transcript, gene_dictionary, downstream_distance, upstream_distance, cage_peaks, dominant=False):
    """
    Score transcript as True or False if there is a CAGE peak within a given distance of the gene.
    Get the chromosome, strand, gene_start from the dictionary
    Turn Gene start into a region start and stop by adding and subtracting the distance.
    Put through the check_for_cage_peak function
    :param transcript:
    :param gene_dictionary:
    :param downstream_distance:
    :param upstream_distance:
    :param cage_peaks:
    :return:
    """
    chromosome, strand, gene_start = gene_dictionary[transcript]
    chromosome = 'chr' + str(chromosome)
    region_start = gene_start - downstream_distance
    region_stop = gene_start + upstream_distance
    if dominant:
        return check_for_dominant_cage_peak(chromosome, strand, region_start, region_stop, cage_peaks)
    else: return check_for_cage_peak(chromosome, strand, region_start, region_stop, cage_peaks)


def describe_pattern_true_and_false(scored_table):
    """

    :param scored_table: Columns ABD and bolean values for each triad
    """
    bool_array = scored_table.to_numpy()
    patterns = np.array([''.join(['T' if cell else 'F' for cell in row]) for row in bool_array])

    # Find unique patterns and their counts
    unique_patterns, counts = np.unique(patterns, return_counts=True)

    # Display the results
    for pattern, count in zip(unique_patterns, counts):
        print(f"{pattern}: {count}")

def list_all_high_confidence_genes_cs():
    gene_list = []
    for chromosome in ['1A', '1B', '1D', '2A', '2B', '2D', '3A', '3B', '3D', '4A',
                       '4B', '4D', '5A', '5B', '5D', '6A', '6B', '6D', '7A', '7B',
                       '7D']:
        genome_annotation = pd.read_csv(f'/Users/mahony/Documents/my_cage_tool_files/CS_annotation_files/'
                                        f'Triticum_aestivum.IWGSC.56.chromosome.{chromosome}.gff3', skiprows=6, sep='\t',
                                        comment='#')

        genome_annotation.columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'ID']
        genome_annotation = genome_annotation[genome_annotation['type'] == 'gene']
        genome_annotation = genome_annotation[
            genome_annotation['chr'] == chromosome]  # I've done this so that fielder and CS are treated the same
        # Strip ID=gene: from the ID column and split at the first ; to remove the rest of the annotation
        genome_annotation['ID'] = genome_annotation['ID'].str.replace('ID=gene:', '')  # ID=gene: for CS
        genome_annotation['ID'] = genome_annotation['ID'].str.split(';').apply(lambda x: x[0])
        gene_list.extend(genome_annotation['ID'].tolist())
    return gene_list


if __name__ == '__main__':

    # gene_dictionary = make_gene_dict(triad_table)
    # with open('/Users/mahony/Documents/1.wheat_triads_project/data/Triad_data/hc_gene_info_dictionary_fielder.pickle',
    #           'wb') as handle:
    #     pickle.dump(gene_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # print("Finished making dict ")
    with open('/Users/mahony/Documents/1.wheat_triads_project/data/Triad_data/hc_gene_info_dictionary.pickle',
              'rb') as handle:
        gene_dictionary = pickle.load(handle)
    # Print some examples of the dictionary


    downstream_distance = 2000
    upstream_distance = 50


    print("######## UNIQUE ONLY #########")
    cage_peaks = pd.read_csv('/Users/mahony/Documents/1.core_wheat/data/nanocage_intermediate_files/post_align_processing_only_unique/para_clustered_tc_2023-07-20_unique_only.csv',
                             sep=',', index_col=0)

    all_hc_genes = list_all_high_confidence_genes_cs()
    gene_dictionary = make_gene_dict_from_list(all_hc_genes)
    print(f"Example of dictionary: {list(gene_dictionary.items())[:5]}")
    all_hc_genes_scored = [score_transcript(transcript, gene_dictionary, downstream_distance, upstream_distance,
                                            cage_peaks) for transcript in all_hc_genes]
    print(f"Number of high confidence genes: {len(all_hc_genes)}")
    print(f"Number of high confidence genes with a CAGE peak: {sum(all_hc_genes_scored)}")
    exit()

    scored_table = triad_table.applymap(lambda transcript: score_transcript(transcript, gene_dictionary,
                                                                       downstream_distance,
                                                                       upstream_distance,
                                                                       cage_peaks))
    print(scored_table.head())
    # Return summary statistics of scored table
    print(scored_table.describe())
    describe_pattern_true_and_false(scored_table)

    print("######## MULTI ONLY #########")
    cage_peaks = pd.read_csv(
        '/Users/mahony/Documents/1.core_wheat/data/nanocage_intermediate_files/post_align_processing_output_multi_unique/para_clustered_tc_2023-07-20_multi_unique.csv',
        sep=',', index_col=0)

    scored_table = triad_table.applymap(lambda transcript: score_transcript(transcript, gene_dictionary,
                                                                            downstream_distance,
                                                                            upstream_distance,
                                                                            cage_peaks))
    print(scored_table.head())
    # Return summary statistics of scored table
    print(scored_table.describe())
    describe_pattern_true_and_false(scored_table)
