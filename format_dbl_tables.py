import csv
import sys
import pandas as pd
import subprocess
from Bio import SeqIO
import re

### Functions

def write_fasta_from_dict(sequence_dict, output_file):
    with open(output_file, 'w') as f:
        for sequence_id, sequence in sequence_dict.items():
            f.write(">{0}\n{1}\n".format(sequence_id, sequence))

def extract_cdhit_sequences_name(input_string):
    unique_sequence_names = re.findall(r'>\w+', input_string)
    unique_sequence_names = [name.strip('>') for name in unique_sequence_names]
    return unique_sequence_names

### Input the sequences from the csv file to a list, update the file with uppercase sequences

Input_table = sys.argv[1]

with open(Input_table, "r") as csv_file:
    table = csv.reader(csv_file)
    next(table)
    seq_column_raw = [row[4] for row in table]

sequences = [string.upper() for string in seq_column_raw]

### Create a dictionary with the varIDs from the list, write fasta file with "unique_sequence" tag

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

reverse_seqdict = {v: k for k, v in seqdict.items()}

### Cluster the sequences and create a dictionary based on their similarities

cmd = ['cd-hit', '-i', 'sequences_unclustered.fasta', '-o', 'output.fasta', '-c', '0.95']

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()

if stderr:
    print(stderr.decode('utf-8'))
else:
    print('CD-HIT completed successfully')

### Import output from CD-hit and parse it

clusters = SeqIO.parse("output.fasta.clstr","fasta")

cluster_lists_ids = []

for record in clusters:
    string = str(record.seq)
    unique_sequences_list = extract_cdhit_sequences_name(string)
    cluster_lists_ids.append(unique_sequences_list)

dbl_dict = {}

for i in range(len(cluster_lists_ids)):
    id_key = "tag-" + str(i+1)
    dbl_dict[id_key] = cluster_lists_ids[i]

reverse_dbl_dict = {}

for k, v in dbl_dict.items():
    for item in v:
        reverse_dbl_dict[item] = k

### Use pandas to update the sequence column and give a unique dbla ID

df = pd.read_csv(Input_table)

df['sequence'] = df['sequence'].str.upper()
df['sequence'] = df['sequence'].map(reverse_seqdict)
df['sequence'] = df['sequence'].map(reverse_dbl_dict)
df = df.rename(columns={'Id':'month','sequence': 'DBL_tag'})

df.to_csv('formatted_table.txt', sep='\t', index=False)