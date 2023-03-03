---
title: README.md
tags: ["Viral metaenomics", "H1N1", "segmented viral genome", "Bioinformatics", "Linux", "Analysis", "Tutorial"]
---
# **Building capacity in Pathogen Genomics: Influenza A Virus - Metagenomics Approach**
---
###### ***Trainers***: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk) & [Gilbert Kibet](https://github.com/kibet-gilbert)
---

- [Introduction](#introduction)
- [Scope of the Tutorial](#scope-of-the-tutorial)
- [Background](#background)
- [Accessing the HPC](#accessing-the-HPC)
- [Bioinformatics Analysis](#Bioinformatics-Analysis)
    - [Step 1 : Preparing the project directory](Step-1---Preparing-the-project-directory)

## Introduction
There are four types of ***influenza viruses*** (A, B, C and D). Influenza A and B cause seasonal flu epidemics in human populations [as reported by CDC](https://www.cdc.gov/flu/about/viruses/types.htm). Influenza C causes mild illness and D infects cattle but has not been detected in humans. Influenza A virus infects birds (avian), humans, pigs (swine), horses (equine), canine (canine) and bats. It is the only one known to cause pandemics (global epidemics), a result of it's ability to mutate quickly and spread efficiently through populations with no or little immunity against it.   

There are four main subtypes of influenza A viruses - ***H1N1***, ***H1N2***, ***H3N2*** and ***H3N1***. The are subtyped based on two surface proteins: ***hemagglutinin (H)*** and ***neuramidase (N)***. Influenza A virus genome is made up of eight segments, H is encoded by segment 4 and N by segment 6. In total there are 18 different hemagglutinin subtypes and 11 different neuraminidase subtypes, with a total of 130 influenza A subtypes so far identified in the wild. Because of this, Influenza virus is prone to reassortment - a process that occure when two influenza subtypes infect the same host concurrently and genetic information (in this case segments) is swapped - creating a new subtype.   
 
Currently two types Influenza (Flu) viruses (A and B), and two influenza A subtypes (H1N1 and H3N2) are routinely detected in different parts of the world. They are classified as shown in the image:
![alt text](https://www.cdc.gov/flu/images/about/influenza-viruses-1200px.jpg?_=72413 "Human Seasonal Influenza Viruses")

## Scope of the Tutorial
We will go through the bioinformatics analysis steps taken when analysing metagenomics sequence data. In a minimal way we will identify the microbes present in the sample. We will then identify our target pathogen and carry out the steps needed to generate the consensus genome, detect the mutations it has and identify the subtype and clade.  

For this tutorial one needs access to a Unix command line in order to access the HPC. One also needs to be proficient in basic linux command line commands.

## Background
In this tutorial we analyse sequences from a clinical sample that were generated in a metagenomics approach, where no targeted amplification of the pathogen was done. The data was generated from a NextSeq 550 at [ILRI](www.ilri.org). The sequencing was done in a `paired-end` approach to give two sequence read files `R1` forward reads and `R2` reverse reads.

## Accessing the HPC
Sign in to HPC using the following command. You will have been assigned a ***username*** that looks like `Bio4InfoXX` and a ***password***.
1. Replace `<user_name>` in the command with the provided username and execute (click enter). 
2. Enter the password and execute. ***Note:*** The password will be typed in the background but will not be visible to you as you type.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
There are two nodes to choose from: `compute05`  and `compute06`. If your username (`Bio4InfoXX`) ends with an ***Odd Number*** (1,3,5,7,9) use `compute05` and if it ends with n ***even number*** (2,4,6,8,0) use `compute06`. Now let us secure a four of CPUs in one of the HPC nodes.  
>Compute05
```
interactive -w compute06 -c 4 -J metagen -p batch
```
>Compute06
```
interactive -w compute06 -c 4 -J metagen -p batch
```

## Bioinformatics Analysis

We will start by seting up the project directory structure and then conduct the analysis stepwise.

### Step 1 : Preparing the project directory
To setup a well-structured project directory we need to create some directories to store our data and scripts. We will be conducting our a anlysis from a directory in the `scratch` space of the HPC.
1. Create a directory using your username in the scratch: 
```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```
2. Create project directories:
```
mkdir -p ilri-africa-cdc-training/viralMetagen/{data,scripts}
cd ilri-africa-cdc-training/viralMetagen/
mkdir -p ./data/{database,fastq,fastqc,fastp,centrifuge,kraken,spades,quast,bowtie,krona,ivar,samtools,snpeff,nextclade}
```

### Step 2: Loading Modules:
Load programs/tools using the following commands:
```
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
module load nextclade/2.11.0
```
Check if modules have been loaded:
```
module list
```
### Step 3: Copying the data and databases:
### Setting up databases
ln -s /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/* ./data/database/
#cp -r /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/* ./data/databaseV1

### Setting up scripts and images
cp -rf /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/scripts/* ./scripts/
1. Copying the data fastq and databases
```
cp /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastq/sample01_R* ./data/fastq/
```

During the analysis we will need databases for different steps. We will explain this during each individual step.
```
cp /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastq/sample01_R* ./data/fastq/
cp -r /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/ ./data/databse
cp -r /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/scripts/* ./scripts/

```

**Note:** *How did we get the fastq files?*

---
<details close>
  <summary>Tip</summary>
  Downloading data from SRA matich the SRR23143759_1<br>
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/059/SRR23143759/SRR23143759_1.fastq.gz -P ./data/fastq<br>
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/059/SRR23143759/SRR23143759_2.fastq.gz -P ./data/fastq<br>
</details>

---
Swine Flu (H1N1) virus infect both pigs, and is the 

### Step 4: Assessing Read Quality using fastqc before quality trimming
```
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/sample01_R1.fastq.gz \
	./data/fastq/sample01_R2.fastq.gz
```

You can proceed and copy fastqc output files to local laptop --- Run this command on your laptop not HPC
```
# scp username@hpc.ilri.cgiar.org:~/ilri-africa-cdc-training/viralMetagen/data/fastqc/*.html ./
```

### Step 5: Quality Trimming fastq files with fastp and Trims adapter sequences
```
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
```

### Step 6: Assessing Read Quality after quality trimming
```
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastp/sample01_R1.trim.fastq.gz \
	./data/fastp/sample01_R2.trim.fastq.gz
```

 Copy the trimmed fastqc output files to local laptop --- Run this command on your laptop not HPC
```
# scp username@hpc.ilri.cgiar.org:~/ilri-africa-cdc-training/viralMetagen/data/fastqc/*.html ./
```

### Step 7 (Optional): Taxonomic Classification of Reads
To quickly profile the taxonomic composition of the reads present in the sequences you can proceed as follows:
```
mkdir data/database/centrifuge/
cd data/database/centrifuge/
```
### Classification
```
centrifuge -x ./data/database/centrifuge/hpvc \
	-1 ./data/fastp/sample01_R1.trim.fastq.gz \
	-2 ./data/fastp/sample01_R2.trim.fastq.gz \
	--report-file ./data/centrifuge/sample01-report.txt \
	-S ./data/centrifuge/sample01-results.txt \
	-p 8 \
	--mm 100GB
```

**Tip:** *How do you Build or access a centrifuge Database?*
---
---

### Visualise the Taxonomic classification results with krona tools
Convert centrifuge report to kraken-like report
```
centrifuge-kreport -x ./data/database/hpvc \
	./data/centrifuge/sample01-results.txt > ./data/centrifuge/sample01-kreport.txt
```
Preparing the classification data
```
cat ./data/centrifuge/sample01-results.txt | cut -f 1,3 > ./data/centrifuge/sample01-results.krona
```
Visiualize the report - create a HTML file
```
apptainer run scripts/singularity/krona_2.7.1--pl526_5.sif \
	ktImportTaxonomy -tax ./data/database/krona/taxonomy \
	-o ./data/centrifuge/sample01-results.html \
	./data/centrifuge/sample01-results.krona > ./data/centrifuge/sample01-results.html
```

**Tip:** *How do you Build or access a centrifuge Database?*
---
---

### Step 8: Filter Host Genome in preparation for genome assembly
```
mkdir ./data/database/host_db
cd ./data/database/host_db
```
### Filtering Host genome seqiuences 
```
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
```
**Tip:** *How do you build or access a host genome database - kraken2*
---
---

## Lets Focus on the target pathogenic virus species: H1N1 - Influenza A Virus

### Step 9: Setting up Reference datasets - reference genome, annotation files
Download Genome from NCBI - Genome database - Reference Genome (Influenza A virus (A/New York/392/2004(H3N2)))
```
mkdir -p ./data/database/refseq/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/085/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna.gz -P ./data/database/refseq/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/085/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff.gz -P ./data/database/refseq/
gunzip data/database/refseq/*.gz
```

Renaming files
```
mv ./data/database/refseq/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna ./data/database/refseq/influenzaA.fna
mv ././data/database/refseq/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff ./data/database/refseq/influenzaA.gff
```

### Step 10: Index reference genome - samtools
```
samtools faidx \
	./data/database/refseq/influenzaA.fna \
	--fai-idx ./data/database/refseq/influenzaA.fna.fai
```

### Step 11: Index reference genome - bowtie
```
mkdir ./data/database/bowtie/
bowtie2-build \
	--threads 4 \
	./data/database/refseq/influenzaA.fna \
	./data/database/bowtie/influenzaA
```

### Step 12: Align reads to reference genome
```
bowtie2 -x ./data/database/bowtie/influenzaA \
	-1 ./data/kraken/sample01.unclassified_1.fastq \
	-2 ./data/kraken/sample01.unclassified_2.fastq \
	--threads 1 \
	--un-conc-gz ./data/bowtie/sample01.unmapped.fastq.gz \
	--local \
	--very-sensitive-local \
	2> ./data/bowtie/sample01.bowtie2.log \
	| samtools view -@ 1 -F4 -bhS -o ./data/bowtie/sample01.trim.dec.bam -
```

### Step 13: Sort and Index aligment map
```
samtools sort -@ 4 \
	-o ./data/bowtie/sample01.sorted.bam \
	-T ./data/bowtie/sample01 \
	./data/bowtie/sample01.trim.dec.bam

samtools index -@ 4 ./data/bowtie/sample01.sorted.bam
```

### Step 13: Coverage computation
```
bedtools genomecov \
	-d \
	-ibam ./data/bowtie/sample01.sorted.bam \
	> ./data/bowtie/sample01.coverage
```

### Step 14: Plot Genome coverage in R
```
Rscript ./scripts/plotGenomecov.R ./data/bowtie/sample01.coverage
```

### Step 15: Consensus Genome construsction
For segmented viruses e.g Influenza A ivar consensus is unable to analyse more than one reference (segment/cromosome) name at once. We need to split by reference:
```
bamtools split -in data/bowtie/sample01.sorted.bam \
	-refPrefix "REF_" \
	-reference
```

Renameing output files
```
rename 'sorted.REF' 'REF' ./data/bowtie/*
```

### Step 16: Loop through segmented BAM files and generate consensus:
```
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
```

### Step 17: Loop through seqmented BAM files and conduct Variant Calling from the alignemnts
```
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
```

### Step 18: Coverting variant files from .tsv to vcf (Variant Call Format) - needed in downstream steps
```
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
```

### Step 19: Annotation of Variants - SnpEff and SnpSift
Annotate the variants VCF file with snpEff
```
for varFile in $(find ./data/ivar/variants -name "*.vcf.gz")
do
	fileName01=`basename -- "$varFile"`
	fileName=${fileName01%.*}
	outName=${fileName%.*}
	java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
		-config ./data/database/snpEff/H1N1/snpEff.config \
		-dataDir ./../ \
		-v H1N1 ${varFile} > ./data/ivar/variants/${outName}.ann.vcf

	# Rename summary.html and genes.txt
	mv ./snpEff_summary.html ./data/ivar/variants/${outName}.ann.summary.html
	mv ./snpEff_genes.txt ./data/ivar/variants//${outName}.ann.genes.txt
	
	#Compress vcf
	bgzip -c ./data/ivar/variants/${outName}.ann.vcf > ./data/ivar/variants/${outName}.ann.vcf.gz
	#Create tabix index - Samtools
	tabix -p vcf -f ./data/ivar/variants/${outName}.ann.vcf.gz
	#Generate VCF files
	bcftools stats ./data/ivar/variants/${outName}.ann.vcf.gz > ./data/ivar/variants/${outName}.ann.stats.txt
done
```
**Tip:** *How do you build a SnpEff Databse?*
---
---

### Step 20: Filter the most significant variants using snpSift
```
for varFile in $(find ./data/ivar/variants -name "*.ann.vcf.gz")
do
	fileName01=`basename -- "$varFile"`
	fileName=${fileName01%.*}
	outName=${fileName%.*}
	java -Xmx4g -jar /export/apps/snpeff/4.1g/SnpSift.jar \
		extractFields \
		-s "," \
		-e "." \
		${varFile} \
		"ANN[*].GENE" "ANN[*].GENEID" \
		"ANN[*].IMPACT" "ANN[*].EFFECT" \
		"ANN[*].FEATURE" "ANN[*].FEATUREID" \
		"ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
		"ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
		"ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
		"ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
		"EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
		> ./data/ivar/variants/${outName}.snpsift.txt
done
```
