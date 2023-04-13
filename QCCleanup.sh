#!/bin/bash

# Input the two paired reads data

Paired_reads1=$1
Paired_reads2=$2
Number_ID1=${Paired_reads1:0:4}
Number_ID2=${Paired_reads2:0:4}

if [ "$Number_ID1" = "$Number_ID2" ] && [ "$Paired_reads1" != "$Paired_reads2" ]; then
  echo "The fastq files are different and have the same prefixes, continuing..."
  echo "ID name: $Number_ID1"
  Number_ID=$Number_ID1
else
  echo "The fastq files have different prefixes or are identical, exiting the shell script..."
  exit
fi

mkdir cleaned_$Number_ID
cp $Paired_reads1 cleaned_$Number_ID/
cp $Paired_reads2 cleaned_$Number_ID/
cd cleaned_$Number_ID

# Run fastQC on both fastq

fastqc $Paired_reads1
fastqc $Paired_reads2

# Remove the first 23 nucleotides and do the using trimmomatic

trimmomatic PE \
 -phred33 $Paired_reads1 $Paired_reads2 fp_$Number_ID.fastq.gz fu_$Number_ID.fastq.gz rp_$Number_ID.fastq.gz ru_$Number_ID.fastq.gz \
 HEADCROP:23 \
 SLIDINGWINDOW:4:20

# Using fastq-join to assemble the paired reads

fastq-join fp_$Number_ID.fastq.gz rp_$Number_ID.fastq.gz -m 4 -p 0 -o joined_$Number_ID.fastq

# Convert fastq files to fasta using seqtk and remove the reads <200bp

seqtk seq -a joined_$Number_ID.fastqjoin > allength_$Number_ID.fa
seqtk seq -L 200 allength_$Number_ID.fa > $Number_ID.fasta

# clean the directory

rm {fu_,ru_,rp_,fp,allength_,joined_}* | rm *{zip,fastq.gz}