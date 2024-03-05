from Bio import SeqIO
import pandas as pd
import sys

def write_fasta_from_dict(sequence_dict, output_file):
    with open(output_file, 'w') as f:
        for sequence_id, sequence in sequence_dict.items():
            f.write(">{0}\n{1}\n".format(sequence_id, sequence))

Sequence_dataframe = sys.argv[1]

allsamples_dataframe = pd.read_csv(Sequence_dataframe)
reduced_dataframe = allsamples_dataframe[['Sample_ID', 'Blast_Cluster_sequnece']]
sequences_dictionary = {}
for index, row in reduced_dataframe.iterrows():
    sequences_dictionary[row['Sample_ID']] = row['Blast_Cluster_sequnece']
write_fasta_from_dict(sequence_dict = sequences_dictionary, output_file = "DBLa_sequences.fa")