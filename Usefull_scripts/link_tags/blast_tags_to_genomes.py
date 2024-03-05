from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import sys
import pandas as pd

#### import files ###

Usage = "Usage: " + sys.argv[0] + "fasta_for_db" + "fasta_to_blast"

if len(sys.argv) != 3:
    print(Usage)
    exit()

Input_database = sys.argv[1]
Input_fasta = sys.argv[2]

### Make the blast database from the first input fasta ###

db_name = "var_genes_database"

makeblastdb_cline = NcbimakeblastdbCommandline(
    dbtype="nucl",  
    input_file=Input_database,
    out=db_name,
)

makeblastdb_cline()

### Set up and launch the blast search

output_file = "blast_var_genes.txt"

blastn_cline = NcbiblastnCommandline(
    query=Input_fasta,
    db=db_name,  
    out=output_file,
    outfmt=6,  # Output format (tabular form)
)

blastn_cline()

#### import the tab-delimited file and extract the sequences with >95% similarity

with open("blast_var_genes.txt", 'r') as tsv_file:
    dictionary_tags = {}
    for line in tsv_file:
        columns = line.split('\t')
        tag = columns[0]
        gene_id = columns[1]
        similarity = float(columns[2])
        if similarity < 95:
            continue
        else:
            if tag in dictionary_tags:                
                dictionary_tags[tag].append(gene_id)
            else:
                dictionary_tags[tag] = [gene_id] ## dont forget to put gene_id as a list

### output the dictionarry in a tidy "tags_varIDs.txt" file

formatted_data = []

for key, values in dictionary_tags.items():
    values_string = str(values)
    values_string = values_string.replace("[","")
    values_string = values_string.replace("]","")
    values_string = values_string.replace("'","")
    formatted_data.append(f"{key},\t{values_string}")

with open("tags_varIDs.txt", 'w') as output:
    for item in formatted_data:
        output.write(item + '\n')

print(dictionary_tags)
