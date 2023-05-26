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

#### format_dbl_tables.py ####

A python script used to compare results of DBLa amplicon sequencing data from multiple batches.

Usage: 

```bash python format_dbl_tables.py Varia_table.csv```

The script will format a csv table with a column containing nucleotide sequences of DBLa clusters, cluster the sequences with a 95% similarity index together (using CDhit) and output a table with newly-generated DBL IDs. See Varia_table.csv for an example of input table, generated from the excel result file of Varia GEM.

#### Varia_excel_parser.R ####

A R script used to parse and syntetise information coming from Varia GEM result output.

Usage: 

```bash Rscript Varia_excel_parser.R```
