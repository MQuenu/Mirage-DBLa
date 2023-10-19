import csv
import sys
import pandas as pd
import subprocess
from Bio import SeqIO
import re

# Define some functions
def write_fasta_from_dict(sequence_dict, output_file):
    # Writes sequences from a dictionary to a FASTA file
    with open(output_file, 'w') as f:
        for sequence_id, sequence in sequence_dict.items():
            f.write(">{0}\n{1}\n".format(sequence_id, sequence))

def extract_cdhit_sequences_name(input_string):
    # Extracts unique sequence names from CD-HIT output
    unique_sequence_names = re.findall(r'>\w+', input_string)
    unique_sequence_names = [name.strip('>') for name in unique_sequence_names]
    return unique_sequence_names

# Input the sequences from the CSV file to a list and convert to uppercase
Input_table = sys.argv[1]

with open(Input_table, "r") as csv_file:
    table = csv.reader(csv_file)
    next(table)  # Skip the header
    seq_column_raw = [row[3] for row in table]

sequences = [string.upper() for string in seq_column_raw]

# Create a dictionary with unique sequence names and write to a FASTA file
count = 0
seqdict = {}

for seq in sequences:
    if seq in seqdict:
        continue
    else:
        count += 1
        VarName = "unique_sequence%i" % (count)
        seqdict[VarName] = seq

write_fasta_from_dict(seqdict, "sequences_unclustered.fasta")

# Reverse the sequence dictionary
reverse_seqdict = {v: k for k, v in seqdict.items()}

# Cluster the sequences using CD-HIT
cmd = ['cd-hit', '-i', 'sequences_unclustered.fasta', '-o', 'output.fasta', '-c', '0.95']

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()

if stderr:
    print(stderr.decode('utf-8'))
else:
    print('CD-HIT completed successfully')

# Import and parse the CD-HIT output
clusters = SeqIO.parse("output.fasta.clstr", "fasta")

cluster_lists_ids = []

for record in clusters:
    string = str(record.seq)
    unique_sequences_list = extract_cdhit_sequences_name(string)
    cluster_lists_ids.append(unique_sequences_list)

# Create a dictionary that maps sequence names to cluster IDs
dbl_dict = {}

for i in range(len(cluster_lists_ids)):
    id_key = "tag-" + str(i+1)
    dbl_dict[id_key] = cluster_lists_ids[i]

# Create a reverse mapping of sequence names to cluster IDs
reverse_dbl_dict = {}

for k, v in dbl_dict.items():
    for item in v:
        reverse_dbl_dict[item] = k

# Use pandas to update the sequence column and rename columns
df = pd.read_csv(Input_table)
df['sequence'] = df['sequence'].str.upper()
df['sequence'] = df['sequence'].map(reverse_seqdict)
df['sequence'] = df['sequence'].map(reverse_dbl_dict)
df = df.rename(columns={'Id': 'month', 'sequence': 'DBL_tag'})

# Write the formatted table to a text file
df.to_csv('formatted_table.txt', sep='\t', index=False)

### Output a fasta file containing the sequences of the sequences linked to a dbla tag

tag_dict = {}

for tags, keys in dbl_dict.items():
    sequence_found = False
    sequences = []
    for key in keys:
        if key in seqdict:
            sequences.append(seqdict[key])
            sequence_found = True
            break
    if sequence_found:
        tag_dict[tags] = sequences

tag_dict_for_output = {key: ' '.join(value) for key, value in tag_dict.items()}

write_fasta_from_dict(tag_dict_for_output, "tags_sequences.fasta")