#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J metagenTest
#SBATCH --qos normal
#SBATCH -n 10

## Changing (switching between) nodes
#interactive -w compute06 -c 10 -J metagen -p batch


## Preparing the project directory:
mkdir -p metagenomics/{data,scripts}
mkdir -p metagenomics/data/{database,fastq,fastqc,fastp,centrifuge,kraken,spades,quast,bowtie,krona,ivar,samtools}
cd metagenomics/

## Downloading data from SRA matich the SRA039136 
## Data from this project: Open-Source Genomic Analysis of Shiga-Toxin–Producing E. coli O104:H4 (https://www.nejm.org/doi/full/10.1056/NEJMoa1107643)
#
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/sample01/sample01_R1.fastq.gz -P ./data/fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/sample01/sample01_R2.fastq.gz -P ./data/fastq

### Checking number of reads in a fastq.gz files:
#zgrep -c '^@' *.fastq.gz
#
### Setting variables for in put
##PROJDIR=$PWD
##FASTQDIR=$PWD/data/fastq/
SAMPLEID=sample01
#
### Loading Modules:
module load fastqc/0.11.9
module load fastp/0.22.0
module load krona/2.8.1
module load centrifuge/1.0.4
module load kraken/2.1.2
module load spades/3.15
module load quast/5.0.2
module load samtools/1.15.1
module load bowtie2/2.5.0
module load bedtools/2.29.0
module load snpeff/4.1g
module load bcftools/1.13
#module load BUSCO/5.2.2
#module load prokka/1.11
#module load bwa/0.7.4
#module load blast/2.12.0+
#
### Assessing Read Quality using fastqc before quality trimming
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/sample01_R1.fastq.gz \
	./data/fastq/sample01_R2.fastq.gz

### Quality Trimming fastq files with fastp and Trims adapter sequences
fastp --in1 ./data/fastq/sample01_R1.fastq.gz \
	--in2 ./data/fastq/sample01_R2.fastq.gz \
	--out1 ./data/fastp/sample01_R1.trim.fastq.gz \
	--out2 ./data/fastp/sample01_R2.trim.fastq.gz \
	--json ./data/fastp/sample01.fastp.json \
	--html ./data/fastp/sample01.fastp.html \
	--failed_out ./data/fastp/sample01_fail.fastq.gz \
	--thread 10 \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--cut_mean_quality 20 \
	--length_required 15 \
	2> ./data/fastp/sample01.fastp.log

### Assessing Read Quality after quality trimming
#
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastp/sample01_R1.trim.fastq.gz \
	./data/fastp/sample01_R2.trim.fastq.gz

## Taxonomic Classification of Reads

# Build Database:
# apptainer pull docker://quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5
# Download NCBI Taxonomy to ./taxonomy/
mkdir data/database/centrifuge/
cd data/database/centrifuge/
centrifuge-download -o taxonomy taxonomy
# Download All complete archaea,bacteria,viral to ./library/
# Downloads from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq availabble domains are: archaea,bacteria,viral,plasmid,fungi,vertebrate_mammalian,vertebrate_other,protozoa,plasmid,plant,metagenomes,mitochondrion,invertebrate,...
centrifuge-download -o library \
	-m \
	-d "archaea,bacteria,viral,plasmid,fungi" refseq > seqid2taxid.map

# Prepare a database - Preffered alternative
wget https://zenodo.org/record/3732127/files/h+p+v+c.tar.gz?download=1
tar -xvzf hpvc.tar.gz 
cd ../../../

#
centrifuge -x ./data/database/centrifuge/hpvc \
	-1 ./data/fastp/sample01_R1.trim.fastq.gz \
	-2 ./data/fastp/sample01_R2.trim.fastq.gz \
	--report-file ./data/centrifuge/sample01-report.txt \
	-S ./data/centrifuge/sample01-results.txt \
	-p 8 \
	--mm 100GB

#Convert centrifuge report to kraken-like report
centrifuge-kreport -x ./data/database/hpvc \
	./data/centrifuge/sample01-results.txt > ./data/centrifuge/sample01-kreport.txt

#Visualization of the taxonomic report using krona
# Load module

#Preparing the data
cat ./data/centrifuge/sample01-results.txt | cut -f 1,3 > ./data/centrifuge/sample01-results.krona
#Build krona db
mkdir ./data/database/krona
apptainer run scripts/singularity/krona_2.7.1--pl526_5.sif \
	ktUpdateTaxonomy.sh ./data/database/krona/taxonomy
apptainer run scripts/singularity/krona_2.7.1--pl526_5.sif \
	ktImportTaxonomy -tax ./data/database/krona/taxonomy \
	-o ./data/centrifuge/sample01-results.html \
	./data/centrifuge/sample01-results.krona > ./data/centrifuge/sample01-results.html

### Filter Host Genome in preparation for genome assembly
## Build host genome database
# Download genome (human)
mkdir ./data/database/host_db
cd ./data/database/host_db
kraken2-build --download-library human \
	--db ./ \
	--threads 4
# Downloading NCBI tax
kraken2-build --download-taxonomy \
	--db ./
# Build database
kraken2-build --build \
	--db ./ \
       	--threads 4
# Removing intermediate files to save space
kraken2-build --clean \
	--db ./
## Alternative - Download prebuilt database
curl -L -o ./kraken2_human_db.tar.gz https://ndownloader.figshare.com/files/23567780
tar -xzvf kraken2_human_db.tar.gz
cd ../../../

#filtering Host genome seqiuences 
kraken2 -db ./data/database/host_db/kraken2_human_db \
	--threads 4 \
	--unclassified-out ./data/kraken/sample01.unclassified#.fastq \
	--classified-out ./data/kraken/sample01.classified#.fastq \
	--report ./data/kraken/sample01.kraken2.report.txt \
	--output ./data/kraken/sample01.kraken2.out \
	--gzip-compressed \
	--report-zero-counts \
	--paired ./data/fastp/sample01_R1.trim.fastq.gz \
	./data/fastp/sample01_R2.trim.fastq.gz

### Genome Assembly using Spades
spades.py -k 27 \
	-1 ./data/kraken/sample01.unclassified_1.fastq \
	-2 ./data/kraken/sample01.unclassified_2.fastq \
	-o ./data/spades/sample01/ \
	-t 8 \
	-m 384
#
###Assess the structure of the genome - examine contiguity
## Run quast
quast.py ./data/spades/sample01/contigs.fasta \
	-t 8 \
	-o ./data/quast/


### Focus on One viral species: H1N1 - Influenza A Virus
# Download Genome from NCBI - Genome database - Reference Genome (Influenza A virus (A/New York/392/2004(H3N2)))
mkdir -p ./data/database/refseq/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/085/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna.gz -P ./data/database/refseq/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/085/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff.gz -P ./data/database/refseq/
gunzip data/database/refseq/*.gz

mv ./data/database/refseq/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna ./data/database/refseq/influenzaA.fna
mv ././data/database/refseq/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff ./data/database/refseq/influenzaA.gff

## Index reference genome - samtools
samtools faidx \
	./data/database/refseq/influenzaA.fna \
	--fai-idx ./data/database/refseq/influenzaA.fna.fai

## Index reference genome - bowtie
mkdir ./data/database/bowtie/
bowtie2-build \
	--threads 4 \
	./data/database/refseq/influenzaA.fna \
	./data/database/bowtie/influenzaA

## Align reads to reference genome
bowtie2 -x ./data/database/bowtie/influenzaA \
	-1 ./data/kraken/sample01.unclassified_1.fastq \
	-2 ./data/kraken/sample01.unclassified_2.fastq \
	--threads 1 \
	--un-conc-gz ./data/bowtie/sample01.unmapped.fastq.gz \
	--local \
	--very-sensitive-local \
	2> ./data/bowtie/sample01.bowtie2.log \
	| samtools view -@ 1 -F4 -bhS -o ./data/bowtie/sample01.trim.dec.bam -

## Sort and Index aligment map
samtools sort -@ 4 \
	-o ./data/bowtie/sample01.sorted.bam \
	-T ./data/bowtie/sample01 \
	./data/bowtie/sample01.trim.dec.bam

samtools index -@ 4 ./data/bowtie/sample01.sorted.bam

## Coverage computation
bedtools genomecov \
	-d \
	-ibam ./data/bowtie/sample01.sorted.bam \
	> ./data/bowtie/sample01.coverage

## Plot Genome coverage in R
Rscript ./scripts/plotGenomecov.R ./data/bowtie/sample01.coverage

## Consensus Genome construsction
# For segmented viruses e.g Influenza A ivar consensus is unable to analyse more than one reference (segment/cromosome) name at once. We need to split by reference:
bamtools split -in data/bowtie/sample01.sorted.bam \
	-refPrefix "REF_" \
	-reference
#Renameing output files
rename 'sorted.REF' 'REF' ./data/bowtie/*

## Loop through segmented BAM files and generate consensus:
mkdir -p ./data/ivar/consensus/
for bamFile in $(find ./data/bowtie -name "*.REF_*.bam")
do
	fileName=`basename -- "$bamFile"`
	outName=${fileName%.*}
	samtools mpileup -aa \
		--count-orphans \
		--no-BAQ \
		--max-depth 0 \
		--min-BQ 0 \
		--reference ./data/database/refseq/influenzaA.fna \
		$bamFile \
		--output ./data/samtools/${outName}.mpileup
	
	cat ./data/samtools/${outName}.mpileup | ivar consensus \
		-t 0.75 \
		-q 20 \
		-m 10 \
		-n N \
		-p ./data/ivar/consensus/${outName}.consensus
done

## Loop through seqmented BAM files and conduct Variant Calling from the alignemnts
mkdir -p ./data/ivar/variants/
for bamFile in $(find ./data/bowtie -name "*.REF_*.bam")
do
	fileName=`basename -- "$bamFile"`
	outName=${fileName%.*}
	samtools mpileup --ignore-overlaps \
		--count-orphans \
		--no-BAQ \
		--max-depth 0 \
		--min-BQ 0 \
		--reference ./data/database/refseq/influenzaA.fna \
		$bamFile \
		--output ./data/samtools/${outName}.var.mpileup
	
	cat ./data/samtools/${outName}.var.mpileup | ivar variants \
		-t 0.25 \
		-q 20 \
		-m 10 \
		-g ./data/database/refseq/influenzaA.gff \
		-r ./data/database/refseq/influenzaA.fna \
		-p ./data/ivar/variants/${outName}.variants
done

## Coverting variant files from .tsv to vcf (Variant Call Format) - needed in downstream steps
for varFile in $(find ./data/ivar/variants -name "*.variants.tsv")
do
	fileName=`basename -- "$varFile"`
	outName=${fileName%.*}
	python3 ./scripts/ivar_variants_to_vcf.py \
		$varFile \
		./data/ivar/variants/${outName}.vcf \
		--pass_only \
		--allele_freq_thresh 0.75 > ./data/ivar/variants/${outName}.counts.log

	#Compress
	bgzip -c ./data/ivar/variants/${outName}.vcf > ./data/ivar/variants/${outName}.vcf.gz
	#Create tabix index - Samtools
	tabix -p vcf -f ./data/ivar/variants/${outName}.vcf.gz
	#Generate VCF files
	bcftools stats ./data/ivar/variants/${outName}.vcf.gz > ./data/ivar/variants/${outName}.stats.txt
done

## Annotation of Variants - SnpEff and SnpSift
for varFile in $(find ./data/ivar/variants -name "*.vcf.gz")
do
	fileName=`basename -- "$varFile"`
	outName=${fileName%.*}
	java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar
done



## Gene content assessment
## Run BUSCO
## Note: to use the full path for output directory, one has to edit the config.ini file and set it from there. Use current dir for now.
#busco \
#-i /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
#-m genome \
#-o /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/busco \
#-l bacteria \
#-c 8
#
## Genome annotation
## Run PROKKA
#prokka \
#/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
#--outdir /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka \
#--cpus 8 \
#--mincontiglen 200 \
#--centre XXX \
#--force
#

## Step 4: Calculate per base coverage with bedtools
#
## First,index the new bam file
#samtools index /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam
#
#bedtools genomecov \
#-ibam /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/sorted_mapped.bam \
#> /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/bwa/coverage.out
#
### AMR identification using RGI: https://github.com/arpcard/rgi#id69
#
##Install/setup RGI- Awaiting instrallation by Alan. For now use singularity image:
## setting up singularity image
#mkdir -p ./scripts/singularity
#apptainer pull ./scripts/singularity/rgi_latest.sif \
#	--force docker://finlaymaguire/rgi:latest
#
## Data Download:
#
#cd ./data/card/
#wget https://card.mcmaster.ca/latest/data -P ./
#wget https://card.mcmaster.ca/latest/variants -P ./
#tar -xvf data ./card.json
#
## Loading Reference data:
#cd ../../
## Check database version - Systemwide
#apptainer run ./scripts/singularity/rgi_latest.sif \
#	rgi database --version
## Clean the loacal database
#apptainer run ./scripts/singularity/rgi_latest.sif \
#	rgi clean --local 
## Load the main database reflecting current version of ref data
#apptainer run ./scripts/singularity/rgi_latest.sif \
#	rgi load --card_json ./data/card/card.json \
#	--local
## Check database version - Local
#apptainer run ./scripts/singularity/rgi_latest.sif \
#	rgi database --version --local
## Perform 
#apptainer run ./scripts/singularity/rgi_latest.sif \
#	rgi main --input_sequence /path/to/nucleotide_input.fasta \
#       	--output_file /path/to/output_file \
#	--local \
#	-a PRODIGAL \
#	--clean \
#	--low_quality \
#	--num_threads 8 \
#	--split_prodigal_jobs
#
