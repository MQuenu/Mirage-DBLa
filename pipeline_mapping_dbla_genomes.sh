#!/bin/bash

input_genome=$1
input_dbla_sequences=$2
prefix_ID="${input_genome:0:7}" 

echo "the genome ID is $prefix_ID"

## Convert the .csv type of files from Varia output to a fasta file

python convert_to_fasta.py $input_dbla_sequences

## map cluster sequences to reference genomes using bwa

bwa index $input_genome
bwa mem $input_genome DBLa_sequences.fa > alignment_dbla_$prefix_ID.sam

## generate separate bam files with aligned / unaligned sequences

samtools view -F 4 alignment_dbla_$prefix_ID.sam > aligned_sequences_$prefix_ID.sam
samtools view -f 4 alignment_dbla_$prefix_ID.sam > unaligned_sequences_$prefix_ID.sam

## Sort and index the sam file

samtools view -bS alignment_dbla_$prefix_ID.sam | samtools sort -o sorted_alignment_dbla_$prefix_ID.bam
samtools index sorted_alignment_dbla_$prefix_ID.bam

## Using mpileup to make a pileup file, then using it to generate .vcf file

samtools mpileup -uf $input_genome -g sorted_alignment_dbla_$prefix_ID.bam > output_$prefix_ID.pileup
bcftools call -mv output_$prefix_ID.pileup > output_$prefix_ID.vcf