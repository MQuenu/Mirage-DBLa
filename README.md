# Mirage-DBLa

This github repository contains an ensemble of scripts that I have created for the analysis of DBLa amplicon sequencing data. It is linked to an 

##### QCCleanup.sh #####

A bash script used to clean up raw DBLa amplicon sequencing data. Uses fastqc v0.12.1, trimmomatics v0.29, fastq-join v1.3.1, seqtk v 1.3.

Usage: 

```bash QCCleanup.sh <forward-reads.fastq> <reverse-reads.fastq>```

Will perform the following steps:

1- Run FastQC on raw reads to assess read quality 

2- Use trimmomatic to trim out first 23bp and poor quality bases 

3- Use fastq-join to join and assemble paired-end reads 

4- convert files from fastq to fasta 

#### Varia_excel_parser.R ####

A R script used to parse the Varia GEM excel spreadsheet into a single tab-delimited file, containing the predicted domains. I found it more useful to work with this than the multi-spreadsheet excel table outputted by Varia.

Usage: 

```Rscript Varia_excel_parser.R <Varia excel results>```

#### tag_sequences.py ####

This python script is used to compare the outputs of different runs of Varia GEM - identify clusters of DBLa common across two different samples.

Usage: 

```python format_dbl_tables.py Varia_table.csv```

The script takes into input a csv table with a column containing the consensus nucleotide sequences of DBLa clusters. This table is typically a concatenation of different varia GEM outputs. The consensus sequences will here be clustered with a 95% similarity index (using CDhit) and output a table where those consensus sequences are replaced by newly-generated DBL tags IDs. See Varia_table.csv for an example of input table. 

Note: the input table can be linked to any metadata that you want, just be sure to have the consensus cluster sequence as row '4' in the input dataframe when running the script


