### Script by Camille Cohen

#!/bin/bash
#Pathway of folders that you use.
GENOME=/Volumes/Camille_Lab/Bio-informatique/MIRAGE_Crick3/01_Reference/plasmodb-46_pfalciparum3d7_genome.fasta
READS=/Volumes/Camille_Lab/Bio-informatique/MGX_Seq_SG/01-Data/
FASTP=/Volumes/Camille_Lab/Bio-informatique/MGX_Seq_SG/04-Fastp/


#Informations avant utilisation
 #1- Ne pas oublier de changer les extensions et nominations des Fastq
 # Versions software
 ## fastqc v0.11.9
 ## fastqscreen v0.14.1
 ## fastp 0.21.0
 ## multiQC 1.9
 ## HISAT2 2.2.1
 ## samtools 1.9
 ## stringtie 2.2.1

#Choose the work volume
cd /Volumes/Camille_Lab/Bio-informatique/MGX_Seq_SG


#FastQC
if [ ! -d "02-FastQC" ];then #Creation du fichier
	mkdir -p 02-FastQC #-p créer tous les dossiers dans l'aborescence
	fastqc -o 02-FastQC -t 4 $READS/*.gz #-o repertoire output pour l'échantillon -t est le nombre de taches que lui dit de faire en mm temps 1 coeur= 1 tache
fi


#FastQScreen
if [ ! -d "03-FastQScreen" ];then #Creation du fichier
	mkdir -p 03-FastQScreen #-p créer tous les dossiers dans l'aborescence
	for ft in $READS/*_R1_001.fastq.gz;do
		r1t=$ft
		r2t=${ft/_R1_001.fastq.gz/}_R2_001.fastq.gz;
		fastq_screen --conf /Volumes/Camille_Lab/Logiciels/FastQ-Screen-0.14.1/fastq_screen_pf.conf --aligner BOWTIE2 $r1t $r2t --outdir 03-FastQScreen #-o repertoire output pour l'échantillon -t est le nombre de taches que lui dit de faire en mm temps 1 coeur= 1 tache
	done
fi


#Fastp

if [ ! -d "04-Fastp" ];then #Creation du fichier
	mkdir -p 04-Fastp #-p créer tous les dossiers dans l'aborescence
	for fq in $READS/*_R1_001.fastq.gz;do
		r1q=$fq
		r2q=${fq/_R1_001.fastq.gz/}_R2_001.fastq.gz;
		echo "Trimming and filtering"
		fastp --in1 $r1q --in2 $r2q --out1 04-Fastp/$(basename $r1q .fastq.gz )_paired.fq.gz --out2 04-Fastp/$(basename $r2q .fastq.gz)_paired.fastq.gz --unpaired1 04-Fastp/$(basename $r1q .fastq.gz)_unpaired.fq.gz --unpaired2 04-Fastp/$(basename $r2q .fastq.gz)_unpaired.fq.gz -w 1 -h 04-Fastp/$(basename $fq)_report_fastp.html -j 04-Fastp/$(basename $fq)_fastp.json
	done
fi


#FastQC after fastp
if [ ! -d "05-FastQC_fastp" ];then #Creation du fichier
	mkdir -p 05-FastQC_fastp #-p créer tous les dossiers dans l'aborescence
	fastqc -o 05-FastQC_fastp -t 4 $FASTP/*.gz #-o repertoire output pour l'échantillon -t est le nombre de taches que lui dit de faire en mm temps 1 coeur= 1 tache
fi


#MultiQC
if [ ! -d "06-MultiQC" ];then #Creation du fichier
	mkdir -p 06-MultiQC #-p créer tous les dossiers dans l'aborescence
	multiQC . -o 06-MultiQC/
	multiQC . -o 03-FastQScreen/*.txt
fi


# Indexing of 3D7 genome (environ 1 min)
if [ ! -d "07-Index-HISAT" ];then #Creation du fichier
	mkdir -p 07-Index-HISAT #-p créer tous les dossiers dans l'aborescence
	hisat2-build $GENOME 07-Index-HISAT/PFalci3D7
fi


# Mapping with HISAT2 on Pf 3D7 (max 2 min par échantillon, total 10 min)
if [ ! -d "08-Mapping-HISAT" ];then #Creation du fichier
	mkdir -p 08-Mapping-HISAT
	for fh in $FASTP/*_R1_001_paired.fq.gz;do
		r1h=$fh
		r2h=${fh/_R1_001_paired.fq.gz/}_R2_001_paired.fastq.gz;
		hisat2 -p 4 --max-intronlen 3000 -x 07-Index-HISAT/PFalci3D7 -1 $r1h -2 $r2h --summary-file 08-Mapping-HISAT/$(basename $fh _R1_001_paired.fq.gz)_summary.txt| samtools view -S -b > ./08-Mapping-HISAT/$(basename $fh _R1_001_paired.fq.gz).bam
	done
fi

#Sorting the BAM files with samtools 

if [ ! -d "09-sortedBAM" ];then #Creation du fichier
	mkdir -p "09-sortedBAM" #-p créer tous les dossiers dans l'aborescence
	for fs in 08-Mapping-HISAT/*.bam;do
		samtools sort $fs > ./09-sortedBAM/$(basename $fs .bam)_sorted.bam
	done
fi

#Filtering the BAM files with samtools 

if [ ! -d "10-FilteredBAM" ];then #Creation du fichier
	mkdir -p "10-FilteredBAM" #-p créer tous les dossiers dans l'aborescence
	for ff in 09-sortedBAM/*.bam;do
		samtools view -b -F 4 $ff > ./10-FilteredBAM/$(basename $ff _sorted.bam)_mapped.bam
		samtools index -M ./10-FilteredBAM/*.bam
	done
fi

#HTseq count. Use a GTF. 


if [ ! -d "11-HTseq" ];then #Creation du fichier
	mkdir -p "11-HTseq" #-p créer tous les dossiers dans l'aborescence
	for fb in ./10-FilteredBAM/*.bam;do
		htseq-count  -t exon -i gene_id -f bam  $fb /Volumes/Camille_Lab/Bio-informatique/GenomeReference/Plasmo3D7/plasmodb-46_pfalciparum3d7_gffread.gtf > ./11-HTseq/$(basename $fb .bam).txt
	done
fi


#MultiQC
if [ ! -d "12-MultiQC" ];then #Creation du fichier
	mkdir -p 12-MultiQC #-p créer tous les dossiers dans l'aborescence
	multiQC . -o 12-MultiQC/

fi