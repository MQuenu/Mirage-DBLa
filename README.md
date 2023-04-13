# Mirage-DBLa

##### QCCleanup.sh #####

A bash script used to clean up raw DBLa amplicon sequencing data. Uses fastqc v0.12.1, trimmomatics v0.29, fastq-join v1.3.1, seqtk v 1.3.

Usage: $bash QCCleanup.sh <forward-reads.fastq> <reverse-reads.fastq>

Will perform the following steps:
1- Run FastQC on raw reads to assess read quality
2- Use trimmomatic to trim out first 23bp and poor quality bases
3- Use fastq-join to join and assemble paired-end reads
4- convert files from fastq to fasta
