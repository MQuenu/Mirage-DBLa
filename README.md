# Mirage-DBLa

##### QCCleanup.sh #####

A bash script used to clean up raw DBLa amplicon sequencing data. Uses fastqc v0.12.1, trimmomatics v0.29, fastq-join v1.3.1, seqtk v 1.3.

Usage: 

```bash QCCleanup.sh <forward-reads.fastq> <reverse-reads.fastq>```

Will perform the following steps:
1- Run FastQC on raw reads to assess read quality
2- Use trimmomatic to trim out first 23bp and poor quality bases
3- Use fastq-join to join and assemble paired-end reads
4- convert files from fastq to fasta

#### tag_sequences.py ####

A python script used to compare results of different Varia GEM runs, and generate DBLa-'tags' that are used to identify shared DBLa sequences across different runs.

Usage: 

```python format_dbl_tables.py Varia_table.csv```

The script will format a csv table with a column containing nucleotide sequences of DBLa clusters, cluster the sequences with a 95% similarity index together (using CDhit) and output a table with newly-generated DBL IDs. See Varia_table.csv for an example of input table, generated from the excel result file of Varia GEM. 

Note: be sure to have the consensus cluster sequence as row '4' in the input dataframe when running the script

#### Varia_excel_parser.R ####

A R script used to parse the Varia GEM result output into a single tab-delimited file

Usage: 

```Rscript Varia_excel_parser.R <Varia excel results>```
