"""
This python script generates a table of fielder triads from a table of CS triads identifying orthologs through
reciprocal blast. Many of these functions I've modifyied from my script cagetool_version2.py

"""

# Import packages
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
import pandas as pd



def find_chromosome_gene_is_on(gene_name):
    """
    Find the chromosome that the gene is on
    :param gene_name: TraesFLD1A01G000200.1 Can be gene or transcript name
    :return: 1A
    """
    # Get the code between the first number and the first G with regex
    # e.g. TraesCS1A02G199600.2 -> 1A
    # e.g. TraesFLD1A01G000200.1 -> 1A
    # e.g. TraesCSU02G202000 -> U
    gene_name = str(gene_name)
    genome_code = re.search(r'Traes(.*?)\d', gene_name).group(1)
    if genome_code == 'CS':
        chromosome = gene_name[7:9]
    elif genome_code == 'FLD':
        chromosome = gene_name[8:10]
    elif genome_code == 'CSU':
        chromosome = gene_name[7]
    return chromosome



def protein_seqrec_cs(transcript_name):
    """
    Get the protein sequence as a SeqRecord object from a CS transcript name.
    :param transcript_name: CS Transcript name e.g. TraesCS1A02G199600.2
    :return: SeqRecord object
    """
    # Path to the FASTA file of all the proteins in the genome

    # At the moment this is hardcoded only for Chinese Spring. In the future this will be changed to having a dictionary
    # of all the genomes and their corresponding protein FASTA files

    chinese_spring_fasta_file = f"input_data/Triticum_aestivum.IWGSC.pep.all.fa"

    # Parse the FASTA file and search for the protein sequence with the corresponding gene name
    for record in SeqIO.parse(chinese_spring_fasta_file, "fasta"):
        if transcript_name in record.description:
            # Create a SeqRecord object with the protein sequence and the corresponding gene information
            protein_seq = Seq(str(record.seq))
            protein_record = SeqRecord(protein_seq, id=record.id, name=record.name, description=record.description)
            return protein_record



def find_fielder_homeolog(gene_name):
    """
    Find the homeolog in Fielder that has the highest blast score
    :param gene_name: TraesCS1A02G199600.2
    :param blast_threshold: 0.9
    :return: TraesFLD1A01G000200.1

    # Note in folder input data I have run the script make_blast_db.sh which performs the following command:
    makeblastdb -dbtype prot -in fielder.release.protein.fa -input_type fasta -blastdb_version 5

    """
    blast_threshold = 0.9
    # Fielder protein FASTA file
    fielder_fasta_file = f"input_data/fielder.release.protein.fa"

    query_gene_protein = protein_seqrec_cs(gene_name)
    # Write this query gene to a FASTA file as having the query as a fasta file seems to work best with
    # NcbiblastpCommandline
    with open(f"intermediate_data/fasta_files_for_blast/{gene_name}.fasta", "w") as output_handle:
        SeqIO.write(query_gene_protein, output_handle, "fasta")

    blastp_cline = NcbiblastpCommandline(query=
                                         f'intermediate_data/fasta_files_for_blast/{gene_name}.fasta',
                                         db=fielder_fasta_file,
                                         outfmt=6,
                                         evalue=blast_threshold)()[0]

    results = blastp_cline.split('\n')
    results = [result.split('\t') for result in results]
    results = pd.DataFrame(results, columns=['Query', 'Subject', '% identity', 'alignment length', 'mismatch',
                                             'gap opens', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    # Remove any empty cells
    results = results.dropna()
    # Filter results for where the ['subject'] column has the same chromosome as the query gene
    input_chromosome = find_chromosome_gene_is_on(gene_name)
    # make a new column in results of the subject chromosome
    results['Subject_chromosome'] = results['Subject'].apply(find_chromosome_gene_is_on)
    # Filter results for where the subject chromosome is the same as the input chromosome
    results = results[results['Subject_chromosome'] == input_chromosome]
    if results.empty:
        return float('NaN')
    else:
        # Sort the results by evalue
        results = results.sort_values(by=['evalue'])
        best_match = results['Subject'][0]
        best_match = re.search(r'(.*?)\.\d', best_match).group(1)
        # if results is an empty dataframe then return NaN
        return best_match

def read_cs_triad_table(hc_triad_file_path):
    """
    Reads in the hc cs triads as a pandas dataframe
    :return: pandas dataframe
    """
    df = pd.read_csv(hc_triad_file_path, sep='\t', header=None, names=['A', 'B', 'D'])
    return df


def generate_fielder_triad_table():
    cs_triad_table = read_cs_triad_table('input_data/hc_triad_table')
    fielder_triad_table = cs_triad_table.applymap(find_fielder_homeolog)
    return fielder_triad_table

if __name__ == '__main__':
    fielder_triad_table = generate_fielder_triad_table()
    fielder_triad_table.to_csv('output_data/fielder_triad_table', sep='\t', index=False, header=False)