---
Title: README.md
Tags: ["Viral metaenomics", "H1N1", "segmented viral genome", "Bioinformatics", "Linux", "Analysis", "Tutorial"]
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
    - [Step 1 : Preparing the project directory](#Step-1--Preparing-the-project-directory)
    - [Step 2: Loading Modules](#Step-2-Loading-Modules)
    - [Step 3: Setting up sequence data and databases](#Step-3-Setting-up-sequence-data-and-databases)
    - [Step 4: Assessing Read Quality using fastqc before quality trimming](#Step-4-Assessing-Read-Quality-using-fastqc-before-quality-trimming)
    - [Step 5: Quality Trimming fastq files with fastp and Trims adapter sequences](#Step-5-Quality-Trimming-fastq-files-with-fastp-and-Trims-adapter-sequences)
    - [Step 6 (Optional): Assessing Read Quality after quality trimming](#Step-6-Optional-Assessing-Read-Quality-after-quality-trimming)
    - [Step 7: Filter Host Genome Reads](#Step-7-Filter-Host-Genome-Reads)
    - [Step 8 (Optional): Taxonomic Classification of Reads](#Step-8-Optional-Taxonomic-Classification-of-Reads)
    - [Visualise the Taxonomic classification results with krona tools](#Visualise-the-Taxonomic-classification-results-with-krona-tools)
- [Lets Focus on a specific pathogenic virus: H1N1 - Influenza A Virus](#Lets-Focus-on-a-specific-pathogenic-virus-H1N1-Influenza-A-Virus)
    - [Step 9: Downloading the Reference data - reference genome and annotation files](#Step-9-Downloading-the-Reference-data-reference-genome-and-annotation files)
    - [Step 10: Indexing the reference genome using samtools and bowtie](#Step-10-Indexing-the-reference-genome-using-samtools-and-bowtie)
    - [Step 11: Align reads to reference genome](#Step-11-Align-reads-to-reference-genome)
    - [Step 12: Sort and Index aligment map](#Step-12-Sort-and-Index-aligment-map)
    - [Step 13: Coverage computation](#Step-13-Coverage-computation)
    - [Step 14: Plot Genome coverage in R](#Step-14-Plot-Genome-coverage-in-R)
    - [Step 15: Consensus Genome construsction](#Step-15-Consensus-Genome-construsction)
    - [Step 16: Loop through segmented BAM files and generate consensus](#Step-16-Loop-through-segmented-BAM-files-and-generate-consensus)
    - [Step 17: Loop through segmented BAM files and call Variants from the alignemnts](#Step-17-Loop-through-segmented-BAM-files-and-call-Variants-from-the-alignemnts)
    - [Step 18: Converting variant files from tsv to vcf (Variant Call Format) - needed in downstream steps](#Step-18-Converting-variant-files-from-tsv-to-vcf-Variant Call Format-needed-in-downstream-steps)
    - [Step 19: Annotation of Variants - SnpEff and SnpSift](#Step-19-Annotation-of-Variants-SnpEff-and-SnpSift)
    - [Step 20: Filter the most significant variants using snpSift](#Step-20-Filter-the-most-significant-variants-using-snpSift)
    - [Step 21: Nextclade Clade assignment](#Step-21-Nextclade-Clade-assignment)

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
interactive -w compute05 -c 4 -J metagen -p batch
```
>Compute06
```
interactive -w compute06 -c 4 -J metagen -p batch
```

## Bioinformatics Analysis

We will start by setting up the project directory structure and then conduct the analysis stepwise.

### Step 1 : Preparing the project directory
To setup a well-structured project directory we need to create some directories to store our data and scripts. We will be conducting our a anlysis from a directory in the `scratch` space of the HPC.
1. *Create a directory using your username in the scratch:*
> **Note:** In this command we use the UNIX environment variable `$USER` which by default was created during logging into the HPC to store your `<user_name>` i.e (`Bio4InfoXX`). You can view its value using the command `echo $USER`.  
```
mkdir -p /var/scratch/$USER
cd /var/scratch/$USER
```
2. *Create project directories:*
> **Note:** We create a project directory `viralMetagen` to store all that pertains to this tutorial/project. Within `viralMetagen` we created `data` and subdirectories to store our input data and results from different analysis steps. We create `scripts` directory to store scripts/code that we genenrate or need in the analysis.
```
mkdir -p ilri-africa-cdc-training/viralMetagen/{data,scripts}
cd ilri-africa-cdc-training/viralMetagen/
mkdir -p ./data/{database,fastq,fastqc,fastp,centrifuge,kraken,spades,quast,bowtie,krona,ivar,samtools,snpeff,nextclade}
```

### Step 2: Loading Modules
Before any analysis we need to load the modules (Bioinformatics programs or tools) we will use in the differemt analyis steps. We can load all modules ahead of the analysis or during every analysis step as we progress with the analysis. We will load all the modules using the following commands:
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
module load bamtools/2.5.1
module load ivar/1.3.1
module load snpeff/4.1g
module load bcftools/1.13
module load nextclade/2.11.0
module load R/4.2
```
> **Note:** `module` is the command used in managing modules. `load` or `list` (see below) are some of it's subcommands. `module load fastqc/0.11.9` will load a tool called `fastqc` and specific version of it `0.11.9` given that we have many versions of `fastqc`.  

Check if modules have been loaded:
```
module list 
```
### Step 3: Setting up sequence data and databases
1. **Copying the data FASTQ:**   
The FASTQ files that we will use for this tutorial have been stored in a directory in the HPC accessible to all:`/var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastq/`. You will need to copy it to the `data/fastq` directory in the `viralMetagen` directory. However, to avoid using so much space you can just create a [`symbolic link`](https://www.futurelearn.com/info/courses/linux-for-bioinformatics/0/steps/201767) to it.
```
ln -s /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastq/sample01_R* ./data/fastq/
```
2. **Setting up the database directories:**   
During the analysis we will need databases to carry out our analysis at different stages. This will be explained in detail during each individual stage where they are needed. For now we will create the directories where we will store the databases.
```
mkdir -p ./data/database{bowtie,centrifuge,hostdb,krona,nextclade,refseq,snpEff}
```

3. **Copying scripts and images:**   
Some of the analysis steps require R and Python scripts. We have them stored in `/var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/scripts`. We will copy this to `./scripts/` within `viralMetagen` as follows:   
```
cp -r /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/scripts/* ./scripts/
```

> **Note:** *You can use data available in NCBI SRA for practice. So how would you get SRA fastq files?*

---
<details close>
  <summary>Tip</summary>
  <blockquote>
    <p dir="auto">
      1. Go to
      <a href="https://www.ncbi.nlm.nih.gov/"
      ref="nofollow">NCBI</a> and search <code>SRA database</code> for 'H1N1'.
    </p>
    <p dir="auto">i
      2. Check the necessary filters: <code>Access: Public</code>, <code>Source: DNA</code>, <code>Type: genome</code>, <code>Library Layout: paired end</code>, <code>Platform: Illumina</code>, <code>File Type: fastq</code>.
    </p>
    <p dir="auto">
      3. Click on a record that was sequenced in a metagenomics workflow
    </p>
    <details close>
      <summary><strong>Example:</strong>Downloading data from SRA: "SRR23143759_1"</summary>
        Forward read: <code>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/059/SRR23143759/SRR23143759_1.fastq.gz -P ./data/fastq</code>
        Reverse read: <code>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/059/SRR23143759/SRR23143759_2.fastq.gz -P ./data/fastq</code>
    </details close>
  </blockquote>
</details>

---

### Step 4: Assessing Read Quality using fastqc before quality trimming
The **raw `FASTQ`** files have sequences as generated by the sequencer. This includes **poor quality reads**, **adapters**, **PhiX** and some reads may be **dublicates**. We need to check the quality of these suquences and clean up if we need to.   
We will use our first module - **FASTQC**, a bioinformatics software used to analyse quality of raw FASTQ format reads and generate visual plots. The report generated by this step will have information on **`Number of reads`,`Sequence Quality`, `Sequence length`, `GC content`, `Adapter content`, `dublication rate`** and others. Run the command below to execute `fastqc` on our reads:

```
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/sample01_R1.fastq.gz \
	./data/fastq/sample01_R2.fastq.gz
```
This will take about 7 Minutes.
Would like to download the results of the `fastqc` command to your local laptop for visualization?.   
If **YES**, you will need to download the `HTML` report to your laptop. Do this by running the following command. 
You can proceed and copy fastqc HTML output files to `~/home`:    
```
mkdir -p ~/viralMetagen
cp /var/scratch/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastqc/*.html ~/viralMetagen/
```
Then run the following command on your local machine.
***Warning!!!:*** ***Run this command on your laptop not HPC***
```
scp <username>@hpc.ilri.cgiar.org:~/viralMetagen/*.html ./
```
> **Discussion:** A report from this step can be found in these links: [sample01_R1](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/fastqc/sample01_R1_fastqc.html) and [sample01_R2](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/fastqc/sample01_R2_fastqc.html). Let us have a discussion about the results as we wait for the run to complete.

### Step 5: Quality Trimming fastq files with fastp and Trims adapter sequences
After assessing the quality, we will proceed and do some Quality Control (QC). With **fastp** module we will trim reads with qualities less than 20 phred score, remove adapters, remove dublicates and remove PhiX if any.
```
fastp --in1 ./data/fastq/sample01_R1.fastq.gz \
	--in2 ./data/fastq/sample01_R2.fastq.gz \
	--out1 ./data/fastp/sample01_R1.trim.fastq.gz \
	--out2 ./data/fastp/sample01_R2.trim.fastq.gz \
	--json ./data/fastp/sample01.fastp.json \
	--html ./data/fastp/sample01.fastp.html \
	--failed_out ./data/fastp/sample01_fail.fastq.gz \
	--thread 4 \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--cut_mean_quality 20 \
	--length_required 15 \
	--dedup \
	|& tee ./data/fastp/sample01.fastp.log
```
This command may take about 15 minutes to complete, so we can proceed and explore the output of the command above in this link: [sample01](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/fastp/sample01.fastp.html)   
You can proceed and copy the fastp HTML output files to local laptop as follows.  
```
cp /var/scratch/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastp/*.html ~/viralMetagen/
```
Then run this command on your laptop not HPC `scp <username>@hpc.ilri.cgiar.org:~/viralMetagen/*.html ./`

### Step 6 (Optional): Assessing Read Quality after quality trimming
After Quality trimming step above we can proceed and check the sequence quality again as described here. Run the command below.
```
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastp/sample01_R1.trim.fastq.gz \
	./data/fastp/sample01_R2.trim.fastq.gz
```
To copy the trimmed fastqc output files to local laptop. Run this commands:
```
cp /var/scratch/gkibet/ilri-africa-cdc-training/viralMetagen/data/fastqc/*trim_fastqc.html ~/viralMetagen/
```
***Run this command on your laptop not HPC***
```
scp <username>@hpc.ilri.cgiar.org:~/viralMetagen/*.html ./
```
> **Discussion:** A report from this step can be found in these links: [sample01_R1.trim](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/fastqc/sample01_R1.trim_fastqc.html) and [sample01_R2.trim](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/fastqc/sample01_R2.trim_fastqc.html). Let us have a discussion about the results as we wait for the run to complete.

### Step 7: Filter Host Genome Reads.
Remember our sequences come from a human host, and it is **NEITHER legal NOR ethical to analyse human sequences without the permision and ethical approval by relevant authorities**. To avoid any such mistakes it is always required that you exclude any human sequences from your reads before any analysis.

In this step we use **`kraken2` module** and a database of human host genome to filter reads mapping to the human genome. *Ideally we are classifing and removing any classified sequences*
> 1. Let us first set up the human host database
```
mkdir -p ./data/database/host_db
ln -s /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/host_db/kraken2_human_db ./data/database/host_db/
```
> 2. Filtering Host genome sequences 
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
> 3. Compress the output of kraken. The fastq files generated by kraken2 need to be compressed to a format (bzip2) useful in the next downstream step.         
> **Note:** This takes aboute 13 minutes
```
bzip2 -fk ./data/kraken/sample01.unclassified_1.fastq ./data/kraken/sample01.unclassified_2.fastq
```
> **Note:** This takes approximately 13 minutes  
#### Building Kraken2 Databases:
##### Alternative 01: Download prebuilt database
```
mkdir -p ./data/database/host_db
cd ./data/database/host_db
curl -L -o ./kraken2_human_db.tar.gz https://ndownloader.figshare.com/files/23567780
tar -xzvf kraken2_human_db.tar.gz
cd ../../../
```
##### Alteranive 02: Build host genome database
> 1. Download genome (human)
```
mkdir -p ./data/database/host_db
cd ./data/database/host_db
kraken2-build --download-library human \
        --db ./ \
        --threads 4
```
> 2. Downloading NCBI tax
```
kraken2-build --download-taxonomy \
        --db ./
```
> 3. Build database
```
kraken2-build --build \
        --db ./ \
        --threads 4
```
> 4. Removing intermediate files to save space
```
kraken2-build --clean \
        --db ./
cd ../../../
```
### Step 8 (Optional): Taxonomic Classification of Reads
Since we are working with metagenomic sequences, it would be fundamental for us
to profile the microbes present in our clinical sample. To do this we need to *`BLAST`*  our sequences against a database of curated genomes whose source organisms are known. There are different databases and accompanying tools for this purpose as can be seen in this publication on [Metagenomic taxonomic classification tools](https://doi.org/10.1016/j.cell.2019.07.010)                            
To quickly profile the taxonomic composition of the reads present in the sequences we chose to work with **[centrifuge](https://ccb.jhu.edu/software/centrifuge/) module**. Centrifuge is a 'rapid and memory-efficient' classification tool for DNA sequences of microbial origin.You can proceed as follows:                 
> 1. set up the reference database:
```
mkdir -p data/database/centrifuge/
cd data/database/centrifuge/
```
> 2. Classification of reads
```
centrifuge -x ./data/database/centrifuge/hpvc \                                
        -1 ./data/kraken/sample01.unclassified_1.fastq.bz2 \                   
        -2 ./data/kraken/sample01.unclassified_2.fastq.bz2 \                   
        --report-file ./data/centrifuge/sample01-report.txt \                  
        -S ./data/centrifuge/sample01-results.txt \                            
        --threads 4 \
        --mm 100GB
```
#### Building centrifuge Database:
##### Alterantive 01: Download prebuilt database
```
mkdir -p data/database/centrifuge/
cd data/database/centrifuge/
wget https://zenodo.org/record/3732127/files/h+p+v+c.tar.gz?download=1
tar -xvzf hpvc.tar.gz
cd ../../../
````
##### Alternative 02: Build a database from NCBI RefSeq  
> 1. Download NCBI `Taxonomy` to `./taxonomy/`
```
mkdir -p data/database/centrifuge/
cd data/database/centrifuge/
centrifuge-download -o taxonomy taxonomy
```
> 2. Download All complete `archaea`,`bacteria`,`viral` and `fungal` to `./library/`  
Downloads can be accessed on `ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq`. Available domains are: `archaea`,`bacteria`,`viral`,`plasmid`,`fungi`,`vertebrate_mammalian`,`vertebrate_other`,`protozoa`,`plant`,`metagenomes`,`mitochondrion`,`invertebrate`,et. cetera
```
centrifuge-download -o library \
        -m \
        -d "archaea,bacteria,viral,plasmid,fungi" refseq > seqid2taxid.map
cd ../../../
```
#### Visualise the Taxonomic classification results with krona tools
> 1. Set up the reference krona Taxonomy database:                             
```
mkdir -p data/database/krona/
ln -s /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/krona/* ./data/database/krona/
```
> 2. Convert centrifuge report to kraken-like report                           
```
centrifuge-kreport -x ./data/database/centrifuge/hpvc \
        ./data/centrifuge/sample01-results.txt > ./data/krona/sample01-kreport.txt
```
> 3. Preparing the classification report data to suitable krona input format.  
> **Note:** *Takes about 40 Seconds*
```
cat ./data/centrifuge/sample01-results.txt | cut -f 1,3 > ./data/krona/sample01-results.krona
```
> 4. Visiualize the report - create a HTML file                                
> **Note:** *Takes about 3 Minutes*
```
apptainer run scripts/singularity/krona_2.7.1--pl526_5.sif \                   
        ktImportTaxonomy -tax ./data/database/krona/taxonomy \                 
        -o ./data/krona/sample01_taxonomy.krona.html \
        ./data/krona/sample01-results.krona 
```
> **Discussion:** A HTML report is generated from this step and can be found in this links: [sample01_taxonomy.krona.html](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/krona/sample01_taxonomy.krona.html)  
#### Building krona database:
Most taxonomic clssification tools use [`NCBI Taxonomy`](https://www.ncbi.nlm.nih.gov/taxonomy). To translate NCBI's taxonomy, a local [taxonomy database](https://github.com/marbl/Krona/wiki/Installing#taxonomy-database) will be needed and can be downloaded from NCBI Taxonomy. 
##### Alterantive 01: Manually Download `taxdump.tar.gz` and build.
NCBI Taxonomy database can be accessed [through`ftp` site](https://ftp.ncbi.nih.gov/pub/taxonomy/) as `taxdump.tar.gz`.
```
mkdir -p ./data/database/krona
cd ./data/database/krona 
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
ktUpdateTaxonomy.sh --only-build ./data/database/krona/taxonomy
cd ../../../
```
##### Alterantive 02: Download build taxonomy database using `updateTaxonomy.sh`
Krona tools has its own scripts to dowload, build or update a taxonomy database from NCBI Taxonomy. The commands to run the scripts are as follows:
```
mkdir -p ./data/database/krona
apptainer run scripts/singularity/krona_2.7.1--pl526_5.sif \
        ktUpdateTaxonomy.sh ./data/database/krona/taxonomy
```
> **Note:** For more information on installing krona databases go to [krona wiki](https://github.com/marbl/Krona/wiki/Installing)

## Lets Focus on a specific pathogenic virus: H1N1 - Influenza A Virus

### Step 9: Setting up Reference datasets - reference genome, annotation files
Based on the taxonomic classification we can tell that we have a high abundance
of ***H1N1 Influenza A Virus***. To conduct analysis focused on this virus we will have to download the refernce Genome data from NCBI Genome database. The Reference Genome for Influenza A virus H1N1 (A/California/07/2009(H1N1))).         
Follow the following steps to identify and download the data:
> 1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
> 2. Type 'H1N1' on the search box and select 'Genome' database. On the landing page you will get you will see:   
  > **Influenza A virus**  
  > **Reference genome:** [Influenza A virus (A/New York/392/2004(H3N2))](https://www.ncbi.nlm.nih.gov/genome/10290?genome_assembly_id=899994)   
 
Because the Influenza A Virus genome segments rapidly evolves, we will select the most recent reference genome to H1N1. This would be A/California/07/2009(H1N1). Reported in the [2009 pandemic](https://doi.org/10.1371/journal.pone.0013381) in California.
> 3. Under [Download sequence and annotation from **RefSeq** or **GenBan*k*] select the [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/).
> 4. Note that there are seven versions of `Influenza_A_virus/latest_assembly_versions`.  
Select the genome version [ViralMultiSegProj274766](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/).  
You can click on `*_assembly_report.txt` to confirm `A/California/07/2009(H1N1)`.
> 5. Right click on the [\*genomic.fna.gz (FASTA)](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.fna.gz) and select 'copy link' and download using `wget` command as follows:
```
mkdir -p ./data/database/refseq/
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.fna.gz -P ./data/database/refseq/
```
> 6. Download the genome annotation GFF file by right clicking on [\*genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.gff.gz) and copying the link then dowloading uisng `wget` command.
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/GCF_001343785.1_ViralMultiSegProj274766_genomic.gff.gz -P ./data/database/refseq/                    
```
> 7. Dowload the [md5checksum](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/md5checksums.txt) and check for integrity of your reference genome (FASTA) and annotation (GFF) files.
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_001343785.1_ViralMultiSegProj274766/md5checksums.txt -P ./data/database/refseq/
cd ./data/database/refseq/
md5sum -c md5checksums.txt | less -S
```
> Change directory back to `viralMetagen`.
```
cd ../../../
```
> 8. Upon successfull integrity check of the files (`OK`), Decompress the `.gz` files
```
gunzip ./data/database/refseq/*.gz
```
> 9. Rename the `FASTA` and `GFF` files
```
rename 'GCF_001343785.1_ViralMultiSegProj274766_genomic' 'H1N1' ./data/database/refseq/*
```

### Step 10: Indexing the reference genome using samtools and bowtie           

> 1. Indexing reference genome `FASTA` using `samtools faidx`.                 
Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname`, `seqlength`, `first-base offset`, `seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for `samtools`' operations and `bcftools` another tool we will later use.                                
```
samtools faidx \
        ./data/database/refseq/H1N1.fna \                                      
        --fai-idx ./data/database/refseq/H1N1.fna.fai                          
```
> **Note:** *Takes less than a second*
> We can view the product of samtools faidx command using the command below. Press `q` to quit.
```
less -s ./data/database/refseq/H1N1.fna.fai                                    
```
> 2. Indexing using `bowtie`.
`bowtie2` will be used in aligning reads to the reference genome. Since our genome has many segments with many nucluotides, there is need to distinctly identify positions in the entire genome by assigning them indices (co-ordinates\*) or rather numbering their positions. To do so run the command below.                
```
mkdir -p ./data/database/bowtie/
bowtie2-build \
        --threads 4 \
        ./data/database/refseq/H1N1.fna \                                      
        ./data/database/bowtie/H1N1
```
> **Note:** *Takes about 2 Seconds*
`bowtie2-build` outputs six files `1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, and `.rev.2.bt2` all constituting the index and are all needed to align reads to the reference which will no longer be needed by bowtie after this indexing. The output is in binary format and can not be visualized like text.        

### Step 11: Align reads to reference genome                                   
As may have been explained earlier, the basis of comparative genomics is *alignment* of reads to the reference genome. Here is where `bowtie2` takes a read and
*lines up* it's characters to a position in the reference genome with the most similar sequence of nucleotides. This is a very difficult computation problem for a number of reasons:
- The reference genome can be big. H1N1 genome is relatively small.
- The reads may not match exactly to its probable locus of origin
- The reads may have some bad quality quality nucleotides
All these and others have to be factored in in the alignment process. [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is ultrafast and memory-efficient but other tools exist e.g. [BWA MEM2](https://github.com/kaist-ina/BWA-MEME)
```
bowtie2 -x ./data/database/bowtie/H1N1 \
        -1 ./data/kraken/sample01.unclassified_1.fastq \
        -2 ./data/kraken/sample01.unclassified_2.fastq \
        --threads 1 \
        --un-conc-gz ./data/bowtie/sample01.unmapped.fastq.gz \
        --local \
        --very-sensitive-local \
        2> ./data/bowtie/sample01.bowtie2.log \
        | samtools view -@ 1 -F4 -bhS -o ./data/bowtie/sample01.trim.dec.bam -
```
> **Note:** *Takes about 35 minutes*

### Step 12: Sort and Index aligment map
The alignment of `FASTQ` files above to the reference genome results in a random arrangement of reads in the order which the sequences occurred in the input FASTQ file. In order to perform any analysis and visualize the alingment, the alignment should be sorted in the order which they occur in the reference genome based on their `alignment coordinates`.
> 1. First let us sort our read alignment `BAM` file above by their occurence in refernce genome.
```
samtools sort -@ 4 \
        -o ./data/bowtie/sample01.sorted.bam \                                 
        -T ./data/bowtie/sample01 \
        ./data/bowtie/sample01.trim.dec.bam                                    
```
> **Note:** Takes about 4 seconds.  
> 2. Then we can now index the genome sorted `BAM` file.  
- Indexing of the BAM file is required by genome viewers like `IGV`, and also by tools that can be used to extract alignments or information like mutation.
- Indexing will generate an `index` file with `.bam.bai` extention.
```
samtools index -@ 4 ./data/bowtie/sample01.sorted.bam
```
> **Note:** *Takes about 5 Seconds*  
- The 'BAM' viewed with `samtools view` has many columns as shown below:
```
NB552490:29:HY2LGBGXK:1:11101:23816:4618        97      NC_026438.1     1      42       16S135M =       39      208     GTCAAATATATTCAATATGGAGAGAATAAAAGAGCTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCATATGGCCATCATCAAAAAGTACACATCGGGAAGGCAAGAGAAGAACCCCGCGCTC AAAAAEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA/EEEE<EEEEEEEEEEEEEEEE<EEEEEEEEEEEEE/EEEEEEAEEAEEEEEA<//EEAEEEEEEE AS:i:242        XN:i:0  XM:i:4  XO:i:0  XG:i:0  NM:i:4 MD:Z:17A71A17A23A3       YS:i:239        YT:Z:DP                                
NB552490:29:HY2LGBGXK:1:11103:24220:6132        97      NC_026438.1     1      42       12S131M1S       =       1       144     AATATATTCAATATGGAGAGAATAAAAGAGCTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCATATGGCCATCATCAAAAAGTACACATCGGGAAGGCAAGAGAAGAACCCCGCG        AAAAAEEEEEEEEEAEEEEEEEEAEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EE/EEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEAEE//AEEAEEE/EEEEEEAAAAEEEEAEEEAEAAEEEEAA</E<EAEEE        AS:i:241        XN:i:0  XM:i:3  XO:i:0  XG:i:0 NM:i:3   MD:Z:17A71A17A23        YS:i:241        YT:Z:DP                        
NB552490:29:HY2LGBGXK:1:11103:18290:12161       161     NC_026438.1     1      42       12S113M =       1       125     AATATATTCAATATGGAGAGAATAAAAGAGCTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCATATGGCCATCATCAAAAAGTACACATCGGGAAG   AAAAAEEEEEEEEEEEEEEEEEEEEE/EAEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEAEEEEEEEEE/EEEEEEEEEAEEAEA<EEEEEEEEEEEEEEEEE<EE   AS:i:205        XN:i:0 XM:i:3   XO:i:0  XG:i:0  NM:i:3  MD:Z:17A71A17A5 YS:i:205        YT:Z:DP        
NB552490:29:HY2LGBGXK:1:11104:4973:14511        97      NC_026438.1     1      42       20S123M =       1       143     GCAGGTCAAATATATTCAATATGGAGAGAATAAAAGAGCTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCATATGGCCATCATCAAAAAGTACACATCGGGAAGGCAAGAGAAG AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEAE AS:i:225        XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:17A71A17A15YS:i:225        YT:Z:DP
NB552490:29:HY2LGBGXK:1:11105:9877:1184 97      NC_026438.1     1       41     16S118M  =       11      134     GTCAAAAATATTCAATATGGAGAGAATAAAAGAGCTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCATATGGCCATCATCAAAAAGTACACATCGGGAAGGCAAG  6AAAAA6EE/AEAEEE/AE/EEEEEEEEEEEEE6EEEEEEEEEE/EE/EEEEEEE/EE/AA6EEAAEEEEAAEEEE//AEEEEA6AEA<6</AEEEEEE//EEEE6E/<AEE/EE/AEEEA<E/AEE/EE/AEE  AS:i:219       XN:i:0   XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:17A71A17A10        YS:i:180       YT:Z:DP
```
This can was earlier explained. The documantation is also available here: [SAM Format](http://samtools.github.io/hts-specs/SAMv1.pdf)                          
> **Quiz:** *To `view` the sorted alignment `BAM` file which command will you use?*
---
<details close>
  <summary>Tip!</summary>
  <blockquote>
    <p dir="auto">
      <code>samtools view</code>
    </p>
    <p dir="auto">
      1. You can view the whole alignment 'BAM' file.
    </p>
    <details close>
      <summary>Answer</summary>
        <code>samtools view ./data/bowtie/sample01.sorted.bam| less -S</code>
    </details close>
    <p dir="auto">
      2. You can count the number of reads in the alignment.
    </p>
    <details close>
      <summary>Answer</summary>
        <code>samtools view ./data/bowtie/sample01.sorted.bam| wc -l</code>
    </details close>
    <p dir="auto">
      3. You can count the number of reads mapping to segment 4 hemagglutinin (HA) gene: NC_026433.1
    </p>
    <details close>
      <summary>Answer</summary>
        <code>samtools view ./data/bowtie/sample01.sorted.bam | grep "NC_026433.1" | wc -l</code>
    </details close>
    <p dir="auto">
      4. You can count the number of reads mapping to segment 6 neuraminidase (NA) gene: NC_026434.1
    </p>
    <details close>
      <summary>Answer</summary>
        <code>samtools view ./data/bowtie/sample01.sorted.bam | grep "NC_026434.1" | wc -l</code>
    </details close>
  </blockquote>
</details>

---

### Step 13: Coverage computation
`bedtools genomecov` computes per-base genome coverage when used with `-d` option. It will do this for all our eight segments.                                 
```
bedtools genomecov \
        -d \
        -ibam ./data/bowtie/sample01.sorted.bam \                              
        > ./data/bowtie/sample01.coverage                                      
```
> **Note:** *Takes about 2 Seconds*  
- The output of this command is a three column file with `chromosome`, `Position` and `depth coverage`. View the ouput with the command below:  
```
less ./data/bowtie/sample01.coverage
```
### Step 14: Plot Genome coverage in R
We can plot the coverage using an R script as shown below. The R script was prepared ahead of the analysis and stored in the `scripts` directory.              
```
Rscript ./scripts/plotGenomeCoverage.R ./data/bowtie/sample01.coverage         
mv ./test_cov_gene_density.png ./data/bowtie/sample01.genomeCoverage.png       
```
The genome coverage of our eight segments can be seen here: [sample01.genomeCoverage](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/bedtools/sample01.genomeCoverage.png)  
> **Note:** *Takes about 5 Seconds*

### Step 15: Consensus Genome construsction                                    
We will use `ivar consensus` to construct the consensus genome.                
Segmented viruses e.g Influenza A virus, are difficult for `ivar consensus` to analyse as it has more than one reference (segment/chromosome). We need to split
the alignment file by reference.
We use `bamtools split` command to split the `BAM` file.                       
```
bamtools split -in ./data/bowtie/sample01.sorted.bam \                         
        -refPrefix "REF_" \
        -reference
```
> **Note:** *Takes about 5 Seconds*
Renaming output files
```
rename 'sorted.REF' 'REF' ./data/bowtie/*                                      
```
Now we have eight bam files for the eight Influenza A virus segments.          
- As we did before, we also need to generate index `.bai` files for these segment `bam` files.
- The command that does this for all `bam` files matching the regex `*REF_*.bam` is shown below:
```
ls ./data/bowtie/*REF_*.bam | xargs -n1 -P5 samtools index -@ 4                
```

### Step 16: Loop through segmented BAM files and generate consensus:
To generate a consensus from the read alignment `BAM` file for one segment `NC_026431.1` we will run the following command:
```
mkdir -p ./data/ivar/consensus/
samtools mpileup -aa \
        --count-orphans \
        --no-BAQ \
        --max-depth 0 \
        --min-BQ 0 \
        --reference ./data/database/refseq/H1N1.fna \
        --region NC_026431.1 \
        ./data/bowtie/sample01.REF_NC_026431.1.bam \
        --output ./data/samtools/sample01.REF_NC_026431.1.mpileup

cat ./data/samtools/sample01.REF_NC_026431.1.mpileup | ivar consensus \
        -t 0.75 \
        -q 20 \
        -m 10 \
        -n N \
        -p ./data/ivar/consensus/sample01.REF_NC_026431.1.consensus
sed -i '/^>/s/Consensus_\(.*\).consensus_threshold.*/\1/' ./data/ivar/consensus/sample01.REF_NC_026431.1.consensus.fa
```
> **Note:** *This takes approximately 2 seconds*
- This has generated a consensus for only one segment: `NC_026431.1`.
- The first part of the command runs `samtools mpileup`to generate a TAB-separated text pileup of the `BAM` files. The format is described in [Pileup Formats](http://www.htslib.org/doc/samtools-mpileup.html#Pileup_Format). Here is a snippet:
```
NC_026431.1     1       A       718     ^M.^M.^M.^M.^K.^M.^K.^M.^M.^M.^K.^K.^M.^M.^K.^M.^K.^M.^M.^M.^M.^K.^M.^M.^K.^M.^M.^M.^M.^M.^M.^K.^K.^M.^K.^M.^M.^M.^E.^M.^M.^M.^K.^M.^K.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^J.^J.^M.^M.^M.^M.^M.^M.^M.^J.^M.^M.^M.^K.^M.^K.^M.^K.^M.^M.^M.^K.^M.^M.^M.^K.^K.^M.^M.^J.^M.^M.^M.^M.^K.^K.^K.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^J.^M.^K.^M.^K.^M.^M.^E.^M.^M.^M.^J.^K.^M.^M.^M.^M.^K.^K.^M.^M.^M.^J.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^K.^M.^M.^K.^K.^M.^M.^M.^M.^E.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^E.^M.^M.^K.^M.^K.^M.^M.^M.^M.^K.^M.^K.^M.^M.^M.^M.^K.^M.^M.^K.^K.^M.^M.^M.^M.^M.^K.^K.^M.^M.^M.^M.^M.^K.^J.^M.^K.^M.^M.^K.^M.^J.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^K.^M.^K.^M.^M.^K.^K.^K.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^K.^K.^K.^M.^K.^M.^M.^M.^K.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^K.^K.^M.^M.^M.^M.^M.^K.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^K.^M.^K.^M.^K.^M.^M.^K.^M.^K.^M.^M.^K.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^K.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^K.^M.^K.^M.^J.^K.^K.^K.^M.^M.^M.^M.^M.^M.^K.^J.^K.^M.^M.^M.^M.^K.^E.^K.^K.^M.^M.^M.^M.^K.^M.^M.^M.^K.^M.^M.^M.^K.^M.^J.^K.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^K.^M.^M.^K.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^E.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^K.^K.^M.^M.^M.^M.^M.^M.^K.^M.^M.^9.^M.^E.^E.^M.^K.^K.^M.^M.^M.^M.^M.^M.^M.^K.^K.^M.^M.^M.^K.^M.^M.^M.^K.^M.^M.^M.^M.^M.^K.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^K,^M,^K,^M,^K,^K,^M,^M,^K,^K,^K,^K,^M,^K,^M,^E,^M,^K,^K,^M,^M,^M,^M,^M,^K,^M,^M,^J,^J,^M,^J,^K,^K,^K,^K,^M,^K,^M,^K,^K,^M,^M,^M,^K,^M,^M,^K,^J,^K,^K,^E,^M,^J,^K,^K,^K,^M,^J,^M,^M,^K,^K,^K,^E,^M,^K,^M,^E,^M,^K,^K,^K,^M,^M,^K,^K,^M,^J,^K,^M,^K,^J,^M,^M,^K,^M,^K,^M,^K,^K,^K,^K,^M,^K,^M,^K,^K,^K,^M,^K,^M,^M,^K,^M,^M,^K,^K,^K,^M,^M,^K,^K,^M,^K,^K,^M,^K,^K,^M,^M,^K,^K,^K,^M,^M,^K,^M,^K,^M,^J,^K,^K,^K,^M,^M,^K,^E,^K,^K,^K,^M,^K,^J,^K,^M,^M,^M,^M,^M,^K,^K,^K,^M,^M,^M,^M,^M,^K,^E,^K,^M,^K,^K,^K,^M,^K,^9,^M,^E,^E,^K,^K,^M,^K,^K,^K,^K,^M,^M,^K,^M,^M,^M,^K,      E/EAEE666EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEAEEEEEEAAEEEAEEEEEEEEEEEE/EEEEEEEEEEEE/EEEEEEAEEEAEEEEEEEEEEEEE/EEEEEEEEEEEAEEEEEEEEEAEEAEEEEEEEEEEEEEEAEEEEEEAA/EAEEE/EEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEAEEEEEE/EE6EEEE/EEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEAEEEEE/EEEEEEEEEEE/EEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEE6AEEEEEEAEEEEEEE6E6EEEEEEEEEEEAEEEEEEEEEEEEEEAEEEEEEEEEE/EE/EEEEEEEEE6EEEAEAEEAEEEAEEAEEE6EEEEEEEEEEEEEE/EEEEEEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEAEAA6/EEAEEEEEEEEE6AEEEEEAA/EEEA/EEAEEE6E6/AAAEEEA/EAEEAEEEEEEAAEEAEE<AEEEE<EEEEEEAAEEEEEEEEEEEE<E<AEEEEAEE/E/EEEEEEEAAEEEEEEAEEAEEEEEEAEEE/EEAEEEEEEAEEEEEEEAE/EAEAEE/EAEEAEE
NC_026431.1     2       T       745     ......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M,^M,^M,^M,^M,^J,^M, EAEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEE/EEEEEEEAEEEAEEAEEEEEEAEEEEEAE/EEAEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEAEEAEEEEEEEAEAEEEEEEEAEEEEEEEEEE/EEEEEEEEEEEEEAEEEEEEEE6/EEEEEE/EAEEEEEEEEEEEEEEEEEEEEEEEAEEEEAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEEEEAE/EEEEEEEEEEEEEEEEEEAE/EE/EEEEEEEAEEEEEEEE/EEEEEEEEEEEEEEEEEEEAEEEEEEEEEEAEEEEAEEEEEEEE/EEEEAE/EEEEEE/EEEEEEEEEE/EEEAEEEEEEEEEEEEEEEEEEEEEAEEE/EEE6EE/EEEEEAEEEEEEEEEAEEEEEEAEEAA/EE/EAEEEEEEE/EAAEEEAEAEEE<EEEAEE<AEAEEEAEEEA/AAEE/EEEEEEEEEAEEEEEEEE<AEEEEEEAEAEEEEEEEEEEEAE<EE<AEEEEAEEEEEEEEEEEEEEAE<<EEEE<EEEAEEEEEAEAEEEAE<EEEEEEEEAE/EEEEEE<EEAE/E/AAAAAAAA!EAAAAA!AAA!EeEaA<e
NC_026431.1     3       G       768     ......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................,,,,,,,^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^M.^K.^M.^M.^M.^M.^M.^M,^M,^M,^M,^K,^M,  EAEAEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEE/E6EEEEEEEEEEEE6EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEE6EEEEEEEEEAAEEEEEEEEEEEEEEEEEE6EAEEEE6EEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEE/EEEEEEEEE6EEEEAEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEAEEEEEE/EEEEEEEEEEEEEEEEE6EAEAEEEEEEEE6EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEEEEAEEAEEEEEEEEEEEEEEEEEEEEEE/EAEEEEEEEEEEEEEEEEEEEEE6EEEEEE/EEEEEEEEEEEEEEEEEEEAEEEEEA/EEEEEEEEEEEEEEEEEEE66EEEE/EEEEEEEEEEEEEEEEEEEAEEEEAEEEEEEEEEEEEEEEEEEEEEE6<EEEAEAEEEEEEEEAEEAEEAAEAEEEEEE<EEE/E</EEEEEEA/EEEE<EEAEEEEEEAEEE<EEEEEAEEEE/E6EEAEE/EEEEEAEEEAEEEEEEEE6EEEEEEEEAEEEEEEEEEEEEEEEEEEEEAEAAE<EEEEE/AEEEEEEEEEAEAEEEEEEEAEAEAAAAAAAA6!EAAAAA!AAA!EeEeEEeAAAAAAAAA!AAAAAAA</<aEA
NC_026431.1     4       A       769     ......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....................,,,,,,,.................,,,,,,^M,     EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/AEEEEEEEEEEEEEEE6E6EEEEEEEEE/EEEEEEEEEEEEEEEEEEEEAEEEEEEEAEEEEEEEEEEEEEEEEEEEEEE/EEEEEEAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEAEEEEEEEEEEAEEAEEEEEEEEEEEEEE/EEEAEEEEEEEEEEEEEEEEEEAEAEEEEEEEE6EEEEEEEEEE/EEEE6E/EEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEEEEEAEEEEEEAEEEEEEAEEEAEEEE6EEEEEEEEEEEEEEEAEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAEEEE/EEEEEEEEEEEEE/E/EEE/EEEEEEEEEEEEEE/EEEEEAEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEEEEEEEEEEEAEEEEAEEEA//EEEAEEEEAAEAEAEEEEEEAAEEEEEEE<AAEE<E6AEEAEEAA/EEEEEEEEEEAEAE<EEEEEEEEE/EAEAEE<EAEEAEEEEEEEEAE/EAEEEAEEAEAEEEEAEEAAEEE/EAEEAEEAEEEEEAAE<EEEEEEAE<EEEEEEEEAEEE/EEEEEAEEEAEEAAAAAAAA!EAAAAA!AAA!EeEeE<eAAAAAAAAA!AA/A/AAAAASE/A
```
- The second part runs `ivar consensus` to convert the `pileup` format to a consensus sequence. The minimum quality accepted is `-q 20`, minimum frequncy of the consensus call is `-t 0.75`, minimum depth is `-m 10` and any position that does not meet the minimum requirements are assinged an `N` i.e `-n N`.            
- The `sed -i` command removes some reduntant phrases from the `FASTA` file header. Changing it from `>Consensus_sample01.REF_NC_026433.1.consensus_threshold_0.75_quality_20` to `>sample01.REF_NC_026433.1`.                                 

- To execute this for all eight segments in a single command, we can run the above command in a loop each time working on a different segment. The following command does this using `for` loops.
```
mkdir -p ./data/ivar/consensus/
for bamFile in $(find ./data/bowtie -name "*.REF_*.bam")
do
        fileName=`basename -- "$bamFile"`
        outName=${fileName%.*}
        chrName=${outName##*REF_}
        samtools mpileup -aa \
                --count-orphans \
                --no-BAQ \
                --max-depth 0 \
                --min-BQ 0 \
                --reference ./data/database/refseq/H1N1.fna \
                --region ${chrName} \
                $bamFile \
                --output ./data/samtools/${outName}.mpileup

        cat ./data/samtools/${outName}.mpileup | ivar consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p ./data/ivar/consensus/${outName}.consensus
        sed -i '/^>/s/Consensus_\(.*\).consensus_threshold.*/\1/' ./data/ivar/consensus/${outName}.consensus.fa
done
```
> **Note:** *Takes about 15 Seconds*
- Consensus sequence files for each of the segments have now been generated and stored in `./data/ivar/consensus/`
- Merge all consensus files for all segments into one file
```
cat ./data/ivar/consensus/*.fa >> ./data/ivar/consensus/sample01.consensus.fa
```

### Step 17: Loop through segmented BAM files and call Variants from the alignemnts
Again we will use the `pileup format` of `BAM` to identify mutations in our samplet. For this we start with segment `NC_026431.1`.
```
mkdir -p ./data/ivar/variants/
cat ./data/samtools/sample01.REF_NC_026431.1.mpileup | ivar variants \
        -t 0.10 \
        -q 20 \
        -m 10 \
        -g ./data/database/refseq/H1N1.gff \
        -r ./data/database/refseq/H1N1.fna \
        -p ./data/ivar/variants/sample01.REF_NC_026431.1.variants
```
> **Note:** *This takes approximately 2 seconds*
- This has generated a file with variants in TSV format for only one segment: `NC_026431.1`.
- The `pileup format` is absolutely the same though with the omission the `-aa` flag in `samtools mpileup` command which only switches off the inclusion of all positions in the reference. See [`samtools mpileup` options](http://www.htslib.org/doc/samtools-mpileup.html#OPTIONS)
- `ivar variants` command calls variants: both single nucleotide variants (SNVs) and indels. `-t 0.10` flag determines the minimum frequency, `-q 20` dictates the minimum quality and `-m 10` the minimum depth.
- The format of the variant call file is TSV as shown below:
```
REGION  POS     REF     ALT     REF_DP  REF_RV  REF_QUAL        ALT_DP  ALT_RV ALT_QUAL ALT_FREQ        TOTAL_DP        PVAL    PASS    GFF_FEATURE     REF_CODON       REF_AA  ALT_CODON       ALT_AA
NC_026431.1     75      G       A       3       1       34      1818    620    34       0.998353        1821    0       TRUE    cds-YP_009118631.1      GCG    AGCA     A
NC_026431.1     75      G       A       3       1       34      1818    620    34       0.998353        1821    0       TRUE    cds-YP_009118628.1      GCG    AGCA     A
NC_026431.1     75      G       A       3       1       34      1818    620    34       0.998353        1821    0       TRUE    cds-YP_009118630.1      GCG    AGCA     A
NC_026431.1     75      G       A       3       1       34      1818    620    34       0.998353        1821    0       TRUE    cds-YP_009121769.1      GCG    AGCA     A
```
- This can be translated [according `ivar` manual](https://andersen-lab.github.io/ivar/html/manualpage.html) as:

|Field | Description |
|----------|:-------------------------------------------|                      
|REGION | Region from BAM file |
|POS | Position on reference sequence |
|REF | Reference base |
|ALT | Alternate Base |
|REF_DP | Ungapped depth of reference base |
|REF_RV | Ungapped depth of reference base on reverse reads |
|REF_QUAL | Mean quality of reference base |
|ALT_DP | Ungapped depth of alternate base |
|ALT_RV | Ungapped deapth of alternate base on reverse reads |
|ALT_QUAL | Mean quality of alternate base |
|ALT_FREQ | Frequency of alternate base |
|TOTAL_DP | Total depth at position |
|PVAL | p-value of fisher's exact test |
|PASS | Result of p-value <= 0.05 |
|GFF_FEATURE | ID of the GFF feature used for the translation |
|REF_CODON | Codong using the reference base |
|REF_AA | Amino acid translated from reference codon |
|ALT_CODON | Codon using the alternate base |
|ALT_AA | Amino acid translated from the alternate codon |

- To execute this for all eight segments in a single command, we will run the command in a loop each time working on a different segment.

```
mkdir -p ./data/ivar/variants/
for bamFile in $(find ./data/bowtie -name "*.REF_*.bam")
do
        fileName=`basename -- "$bamFile"`
        outName=${fileName%.*}
        chrName=${outName##*REF_}
        cat ./data/samtools/${outName}.mpileup | ivar variants \
                -t 0.25 \
                -q 20 \
                -m 10 \
                -g ./data/database/refseq/H1N1.gff \
                -r ./data/database/refseq/H1N1.fna \
                -p ./data/ivar/variants/${outName}.variants
done
```
> **Note:** *This takes approximately 12 seconds*

Variant call files for all segments have been generated and stored in `./data/ivar/consensus/`


### Step 18: Converting variant files from tsv to vcf (Variant Call Format) - needed in downstream steps
The standard used often for reporting and analysing variant calls is [**`Variant Call Format`**](https://samtools.github.io/hts-specs/VCFv4.1.pdf) in short **`VCF`**
We will therefore convert our variant `TSV` to `VCF`. Below we use a python script to perform this: `python3 ./scripts/ivar_variants_to_vcf.py`
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
- You will note that we also compressed the resulting `VCF` file to `.vcf.gz` using `bgzip -c`. This is to prepare it for follow-up analysis.
- Also we index the compressed `.vcf.gz` using `tabix -p` command.
- `bcftools stats` generates some summary statistics from each indexed `vcf.gz` that be used later in `MALTIQC` analysis.

> **Note:** *Takes about 4 Seconds*


### Step 19: Annotation of Variants - SnpEff and SnpSift
1. Set up the SnpEff database first:
```
mkdir -p ./data/database/snpEff/
cp -rf /var/scratch/global/gkibet/ilri-africa-cdc-training/viralMetagen/data/database/snpEff/ ./data/database/snpEff/
```
2. We can now annotate the variants in the compressed `vcf.gz` file with snpEff

```
for varFile in $(find ./data/ivar/variants -name "*.vcf.gz")
do
        fileName01=`basename -- "$varFile"`
        fileName=${fileName01%.*}
        outName=${fileName%.*}
        #SnpEff - annotation
        java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
                -config ./data/database/snpEff/H1N1/snpEff.config \
                -dataDir ./../ \
                -v H1N1 ${varFile} > ./data/ivar/variants/${outName}.snpeff.vcf

        # Rename summary.html and genes.txt
        mv ./snpEff_summary.html ./data/ivar/variants/${outName}.snpeff.summary.html
        mv ./snpEff_genes.txt ./data/ivar/variants//${outName}.snpeff.genes.txt

        #Compress vcf
        bgzip -c ./data/ivar/variants/${outName}.snpeff.vcf > ./data/ivar/variants/${outName}.snpeff.vcf.gz
        #Create tabix index - Samtools
        tabix -p vcf -f ./data/ivar/variants/${outName}.snpeff.vcf.gz
        #Generate VCF files
        bcftools stats ./data/ivar/variants/${outName}.snpeff.vcf.gz > ./data/ivar/variants/${outName}.snpeff.stats.txt
done
```
- `snpEff` also generates summary HTML files stored in `./data/ivar/variants/sample01.REF_NC_*.variants.ann.summary.html`.                                     
- Here are the summary HTML files for all [the eight segments](https://hpc.ilri.cgiar.org/~gkibet/ilri-africa-cdc-training/snpEff/).                           

#### Building SnpEff Database:
##### Alternative 01 - Build Database from `GFF` and `FNA`  
- Set up the Genome Feature File `GFF` and Reference Genome `fna`  
- Create a [`snpEff.config`](https://pcingola.github.io/SnpEff/examples/#step-1-build-database) file in the [right format](https://pcingola.github.io/SnpEff/se_build_db/#add-a-genome-to-the-configuration-file).
```
mkdir -p ./data/database/snpEff/H1N1/
cp ./data/database/refseq/H1N1.gff ./data/database/snpEff/H1N1/genes.gff
cp ./data/database/refseq/H1N1.fna ./data/database/snpEff/H1N1/sequences.fa
echo -e "# Influenza A virus genome, version influezaA\nH1N1.genome: H1N1" > ./data/database/snpEff/H1N1/snpEff.config
```
- Now build a SnpEff database:
```
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./data/database/snpEff/H1N1/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v H1N1
```
##### Alternative 02 - Download Pre-built database:
- Check all pre-built databases (over 20,000) for pathogen of interest.
```
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar databases > ./data/database/snpEff/snpeff.databases.txt
```
View the list of available databases:
```
less -S ./data/database/snpEff/snpeff.databases.txt
```
- Get the genome version and download target
```
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar download -v <genome_version>
```

### Step 20: Filter the most significant variants using snpSift   
[**`SnpSift`**](https://pcingola.github.io/SnpEff/ss_introduction/) is a toolbox with many utilities: `filter`,`annotate`, `caseControl`, `extractFields`, `intervals` and a few others.
- This tutorial uses `extractFields` operation to 'Extract fields from a VCF file to a TXT (tab separated) format'.
- Run the following command to get the most relevant fields.
```
for varFile in $(find ./data/ivar/variants -name "*.snpeff.vcf.gz")            
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
> **Note:** *Takes about 6 Seconds*
- The output of `snpSiff extractfields` command looks like the snippet below:
```
#ANN[*].GENE    ANN[*].GENEID   ANN[*].IMPACT   ANN[*].EFFECT   ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE  ANN[*].RANK     ANN[*].HGVS_C   ANN[*].HGVS_P  ANN[*].CDNA_POS  ANN[*].CDNA_LEN ANN[*].CDS_POS  ANN[*].CDS_LEN  ANN[*].AA_POS  ANN[*].AA_LEN    ANN[*].DISTANCE EFF[*].EFFECT   EFF[*].FUNCLASS EFF[*].CODON   EFF[*].AA        EFF[*].AA_LEN
M1,M2   gene-UJ99_s7gp2,gene-UJ99_s7gp1 LOW,MODIFIER    synonymous_variant,intron_variant       transcript,transcript   Transcript_gene-UJ99_s7gp2,Transcript_gene-UJ99_s7gp1   Coding,Coding   1,1     c.75G>A,c.26+49G>A      p.Ala25Ala,.   75,-1    759,-1  75,-1   759,-1  25,-1   252,-1  0,0     synonymous_variant,intron_variant       NONE,NONE       c.75G>A,c.26+49G>A      p.Ala25Ala,.    252,-1 
M1,M2   gene-UJ99_s7gp2,gene-UJ99_s7gp1 LOW,MODIFIER    synonymous_variant,intron_variant       transcript,transcript   Transcript_gene-UJ99_s7gp2,Transcript_gene-UJ99_s7gp1   Coding,Coding   1,1     c.102A>G,c.26+76A>G     p.Gly34Gly,.   102,-1   759,-1  102,-1  759,-1  34,-1   252,-1  0,0     synonymous_variant,intron_variant       NONE,NONE       c.102A>G,c.26+76A>G     p.Gly34Gly,.    252,-1
M1,M2   gene-UJ99_s7gp2,gene-UJ99_s7gp1 MODERATE,MODIFIER       missense_variant,intron_variant transcript,transcript   Transcript_gene-UJ99_s7gp2,Transcript_gene-UJ99_s7gp1   Coding,Coding   1,1     c.124C>A,c.26+98C>A     p.Leu42Ile,.   124,-1   759,-1  124,-1  759,-1  42,-1   252,-1  0,0     missense_variant,intron_variant NONE,NONE       c.124C>A,c.26+98C>A     p.Leu42Ile,.    252,-1
M1,M2   gene-UJ99_s7gp2,gene-UJ99_s7gp1 MODERATE,MODIFIER       missense_variant,intron_variant transcript,transcript   Transcript_gene-UJ99_s7gp2,Transcript_gene-UJ99_s7gp1   Coding,Coding   1,1     c.238G>A,c.26+212G>A    p.Val80Ile,.   238,-1   759,-1  238,-1  759,-1  80,-1   252,-1  0,0     missense_variant,intron_variant NONE,NONE       c.238G>A,c.26+212G>A    p.Val80Ile,.    252,-1
```
- The result of `snpSift extractFields` can be best understood through this [Documantation page](https://pcingola.github.io/SnpEff/ss_extractfields/).

### Step 21: Nextclade Clade assignment
[**Nextclade**](https://docs.nextstrain.org/projects/nextclade/en/stable/) is a
tool within the [**Nextrain**](https://nextstrain.org/) collection that performs pathogen centered phylogeograhic analysis.                                    
- It uses sequence differences in Multiple Sequence Alignment to assign to [clades](https://clades.nextstrain.org/).
- It also reports suspect quality issues with such sequences.                  
- There are both [web-](https://clades.nextstrain.org/) and [command-line-interfaces](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) for *nextclade*.
To run it in the command-line, we need some reference files: `genome`, `feature map`, `origin tree`, `primers` if used and `quality configurations`. Luckily, for H1N1 and a few other viral pathogens, these files have already been created and can be easily retrieved using the same tool. Otherwise, you will have to create/retrieve accordingly.

> 1. Get the reference data

- Check the list of available nextclade datasets:
```
mkdir -p ./data/database/nextclade/H1N1/
nextclade dataset list > ./data/database/nextclade/nextcladelist.txt
```
- Identify the relevant H1N1 strain based on the reference genome we used.
```
less -S ./data/database/nextclade/nextcladelist.txt
```
- Select and Download latest H1N1 datasets
```
nextclade dataset get \
        --name 'flu_h1n1pdm_ha' \
        --reference 'CY121680' \
        --output-dir ./data/database/nextclade/H1N1/
```
> **Note** *This command takes approximately 10 seconds to complete*           
- You have now downloaded Nextclade reference data for this analysis.          
- View the dataset using `ls ./data/database/nextclade/`                       
- Nextclade datasets are updated regularly and it is advisable to rerun tje download command before any analysis. This will ensure usage of up-to-date datasets.

> 2. Perform Clade assignment:

Now execute nextclade analysis on your sequences.                              
- If you paid attention to the `./data/database/nextclade/nextcladelist.txt` tje you will remamber seeing a dataset name like `flu_h1n1pdm_ha`. The `ha` in the
name stands for `hemagglutinin (H)`. This is because Nextclade analysis only relys on the `hemagglutinin (H)` segment alone. That is `REF_NC_026433.1`.        
- Proceed with your analysis providing the consensus matching `REF_NC_026433.1`.```
nextclade run \
        data/ivar/consensus/sample01.REF_NC_026433.1.consensus.fa \            
        --input-dataset ./data/database/nextclade/H1N1/ \                      
        --output-tsv ./data/nextclade/sample01.REF_NC_026433.1.nextclade.tsv \ 
        --output-tree ./data/nextclade/sample01.REF_NC_026433.1.nextclade.auspice.json \
        --output-all ./data/nextclade/ \                                       
        --output-basename sample01.REF_NC_026433.1.nextclade \                 
        |& tee ./data/nextclade/sample01.REF_NC_026433.1.nextclade.log         
```
- The outcome of this analysis will be stored in `./data/nextclade/`           
> 3. Visualization: The output of Nextclade includes a phylogenetic tree in `.json` format. This tree can be visualized in [Auspice](https://auspice.us/).     

- First let us download the the `.json` file from the HPC:                     
```
cp ./data/nextclade/*.nextclade.auspice.json ~/                                
```
- In your local computer use `scp` to copy the file to any desired destination:
```
scp <user_name>@hpc.ilri.cgiar.org:~/*.nextclade.auspice.json <destination_folder>
```
- Open [Auspice](https://auspice.us/).
- Drag and drop the `.json` file in the [Auspice](https://auspice.us/).        
- Now edit the tree.
  - In `Dataset` click the drop down arrow and select `ncov`, below it select `open` and below it select `global`.
  - In `Color By` click the drop down arrow and select `clade`.                
  - Do any other adjustments as you wish.                                      
                                                  
