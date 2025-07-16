#seqnames	start	end	width	strand	annotation	genes
#Chr1A	211	211	1	+	unknown	
#Chr1A	265	265	1	+	unknown	
#Chr1A	465	465	1	+	unknown	
#Chr1A	1017	1017	1	+	unknown	
#Chr1A	1246	1246	1	+	unknown	
#Chr1A	1340	1340	1	+	unknown	
#Chr1A	5168	5168	1	+	unknown	
#Chr1A	5236	5236	1	+	unknown	
#Chr1A	5634	5634	1	+	unknown	


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Read the text file into a pandas DataFrame
df = pd.read_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/up_500_down_500_shared_R_cadenza_annotated_hc_genes_ctss.csv', sep='\t', comment='#')

# Filter rows where annotation is 'promoter'
promoters = df[df['annotation'] == 'promoter']

# Function to extract sequence
def extract_sequence(row, genome_seq):
    start = row['start'] - 150
    end = row['start'] + 10
    sequence = genome_seq[row['seqnames']].seq[start:end]
    if row['strand'] == '-':
        sequence = sequence.reverse_complement()
    return SeqRecord(sequence, id=f"{row['seqnames']}:{start}-{end}({row['strand']})", description="")

# Load the genome sequence (assuming it's in a FASTA file)
# Genome data  - genome.fasta
genome_fasta = '/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_assembly.fa'
genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))

# Extract sequences
sequences = [extract_sequence(row, genome_seq) for index, row in promoters.iterrows()]

# Write to a new FASTA file
SeqIO.write(sequences, 'promoters.fasta', 'fasta')