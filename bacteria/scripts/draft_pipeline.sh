#!/usr/bin/env bash


## Changing (switching between) nodes
interactive -w compute07 -c 8 -J bact-Test -p highmem


## Preparing the project directory:
mkdir -p bacteria/{data,scripts}
mkdir -p bacteria/data/{fastq,fastqc,fastp,spades}
cd bacteria/

# Downloading data from SRA matich the SRA039136 
# Data from this project: Open-Source Genomic Analysis of Shiga-Toxin–Producing E. coli O104:H4 (https://www.nejm.org/doi/full/10.1056/NEJMoa1107643)

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_1.fastq.gz -P ./data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_2.fastq.gz -P ./data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292862/SRR292862_1.fastq.gz -P ./data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292862/SRR292862_2.fastq.gz -P ./data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292678/SRR292678_1.fastq.gz -P ./data/fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292678/SRR292678_2.fastq.gz -P ./data/fastq

## Checking number of reads in a fastq.gz files:
zgrep -c '^@' *.fastq.gz

## Setting variables for in put
#PROJDIR=$PWD
#FASTQDIR=$PWD/data/fastq/
#SAMPLEID=SRR292770

## Loading Modules:
module load fastqc/0.11.9
module load fastp/0.22.0
module load spades/3.15
module load quast/5.0.2
module load BUSCO/5.2.2
module load prokka/1.11


## Assessing Read Quality using fastqc before quality trimming
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/SRR292770_1.fastq.gz \
	./data/fastq/SRR292770_2.fastq.gz

## Quality Trimming fastq files with fastp
fastp --in1 ./data/fastq/SRR292770_1.fastq.gz \
	--in2 ./data/fastq/SRR292770_2.fastq.gz \
	--out1 ./data/fastp/SRR292770_1.trim.fastq.gz \
	--out2 ./data/fastp/SRR292770_2.trim.fastq.gz \
	--json ./data/fastp/SRR292770.fastp.json \
	--html ./data/fastp/SRR292770.fastp.html \
	--failed_out ./data/fastp/SRR292770_fail.fastq.gz \
	--thread 10 \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--cut_mean_quality 20 \
	--length_required 15 \
	2> ./data/fastp/SRR292770.fastp.log

## Assessing Read Quality after quality trimming

fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastp/SRR292770_1.trim.fastq.gz \
	./data/fastp/SRR292770_2.trim.fastq.gz

## Genome Assembly using Spades
spades.py -k 27 \
	-1 ./data/fastp/SRR292770_1.trim.fastq.gz \
	-2 ./data/fastp/SRR292770_2.trim.fastq.gz \
	-o ./data/spades/ \
	-t 8 \
	-m 384

##Assess the structure of the genome - examine contiguity
# Run quast
quast.py \
/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
-t 8 \
-o /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/quast


# Gene content assessment
# Run BUSCO
# Note: to use the full path for output directory, one has to edit the config.ini file and set it from there. Use current dir for now.
busco \
-i /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
-m genome \
-o busco-results \
-l bacteria \
-c 8

# Genome annotation
# Run PROKKA
prokka \
/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
--outdir /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka \
--cpus 8 \
--mincontiglen 200 \
--centre XXX \
--force

# PROKKA had errors while testing. Will use Glimmer instead

