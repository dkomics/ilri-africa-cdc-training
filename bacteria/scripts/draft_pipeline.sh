#!/usr/bin/env bash


## Changing (switching between) nodes
interactive -w compute07 -c 8 -J bact-Test -p highmem


## Preparing the project directory:
mkdir -p bacteria/{data,scripts}
mkdir -p bacteria/data/{fastq,fastqc,fastp,spades,card}
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
module load bwa/0.7.4
module load samtools/1.9
module load bedtools/2.29.0
module load blast/2.12.0+

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
-o /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/busco \
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

# Extract all the annotated products from the GFF
grep -o "product=.*" /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka/PROKKA_*.gff \
| sed 's/product=//g' | sort | uniq -c | sort -nr > \
/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka/protein_abundances.txt

# Extract the 16S gene sequence
# NB: Can be used for organism identification
./extract_sequences.sh "16S ribosomal RNA" \
/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka/*.ffn \
> /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka/16S_sequence.fasta


# Organism identification
# Using blast to search against the entire nucleotide database.
./blob_blast.sh /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta


# Read mapping to determine coverage
fasta=/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/ecoli_genome/Escherichia_coli_str._K-12_genome.fasta
forward=/var/scratch/global/gkibet/ilri-africa-cdc-training/bacteria/data/fastp/SRR292770_1.trim.fastq.gz
reverse=/var/scratch/global/gkibet/ilri-africa-cdc-training/bacteria/data/fastp/SRR292770_1.trim.fastq.gz


# Step 1: Index your reference genome. This is a requirement before read mapping.
bwa index $fasta

# Step 2: Map the reads and construct a SAM file.
bwa mem -t 8 $fasta $forward $reverse > /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/raw_mapped.sam

# Step3: Remove sequencing reads that did not match to the assembly and convert the SAM to a BAM.
samtools view -@ 8 -Sb  /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/raw_mapped.sam \
| samtools sort -@ 8 -o /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam

# Examine how many reads mapped with samtools
samtools flagstat /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam

# Step 4: Calculate per base coverage with bedtools

# First,index the new bam file
samtools index /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam

bedtools genomecov \
-ibam /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam \
> /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/coverage.out

## AMR identification using RGI: https://github.com/arpcard/rgi#id69

#Install/setup RGI- Awaiting instrallation by Alan. For now use singularity image:
# setting up singularity image
mkdir -p ./scripts/singularity
apptainer pull ./scripts/singularity/rgi_latest.sif \
	--force docker://finlaymaguire/rgi:latest

# Data Download:

cd ./data/card/
wget https://card.mcmaster.ca/latest/data -P ./
wget https://card.mcmaster.ca/latest/variants -P ./
tar -xvf data ./card.json

# Loading Reference data:
cd ../../
# Check database version - Systemwide
apptainer run ./scripts/singularity/rgi_latest.sif \
	rgi database --version
# Clean the loacal database
apptainer run ./scripts/singularity/rgi_latest.sif \
	rgi clean --local 
# Load the main database reflecting current version of ref data
apptainer run ./scripts/singularity/rgi_latest.sif \
	rgi load --card_json ./data/card/card.json \
	--local
# Check database version - Local
apptainer run ./scripts/singularity/rgi_latest.sif \
	rgi database --version --local
# Perform 
apptainer run ./scripts/singularity/rgi_latest.sif \
	rgi main --input_sequence /path/to/nucleotide_input.fasta \
       	--output_file /path/to/output_file \
	--local \
	-a PRODIGAL \
	--clean \
	--low_quality \
	--num_threads 8 \
	--split_prodigal_jobs

