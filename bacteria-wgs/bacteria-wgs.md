---
title: README.md
tags: [ "Pathogen Genomics", "Bioinformatics", "Metadata", "Linux", "Analysis", "Activity"]
---
# **Building capacity in pathogen genomics in Africa**
---
###### ***Trainers***: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk) & [Gilbert Kibet](https://github.com/kibet-gilbert)
---

- [**Building capacity in pathogen genomics in Africa**](#building-capacity-in-pathogen-genomics-in-africa)
          - [***Trainers***: John Juma, Kennedy Mwangi \& Gilbert Kibet](#trainers-john-juma-kennedy-mwangi--gilbert-kibet)
  - [Introduction](#introduction)
  - [Background](#background)
  - [Scope](#scope)
  - [Prerequisite](#prerequisite)
    - [Set-Up](#set-up)
    - [Preparations](#preparations)
      - [***Log into the HPC***](#log-into-the-hpc)
      - [***Project organisation***](#project-organisation)
      - [***Data retrieval and integrity checks***](#data-retrieval-and-integrity-checks)
  - [Analysis](#analysis)
      - [***Loading modules***](#loading-modules)
      - [***Prepare the reference genome***](#prepare-the-reference-genome)
      - [***Quality assessment***](#quality-assessment)
      - [***Quality and adapter filtering***](#quality-and-adapter-filtering)
      - [***Perform genome assembly***](#perform-genome-assembly)
      - [***Genome assessment***](#genome-assessment)
      - [***Genome annotation***](#genome-annotation)
      - [***Organism identification***](#organism-identification)
      - [***Identification of virulence factors***](#identification-of-virulence-factors)
      - [***AMR identification***](#amr-identification)


## Introduction
*E. coli* are ubiquitous bacteria found in the environment including the gut of humans and other animals and consists of numerous types and strains. Most types of *E. coli* are normal inhabitants of the gut of animals and do not cause diseases. However, some of the strains of *E. coli* have acquired genes that enable them to cause disease. These strains are commonly associated with food poisoning leading to diarrhoea and are referred to as diarrheagenic *E. coli* (DEC). Transmission occurs primarily through contaminated food but can also occur via person-to-person transmission, animal contact and water. Shiga toxin-producing *E. coli* (STEC) serotype O157:H7 causes bloody diarrhoea and has previosuly been responsible for outbreaks worldwide.  

> **Note**

> This training workshop was organized by the Africa Centers for Disease Control & Prevention (Africa CDC) jointly with the African Society for Laboratory Medicine (ASLM) and the International Livestock Research Institute (ILRI) in Kenya to build capacity in pathogen genomic surveillance in Africa.


## Background
In May and June of 2011, an *E. coli* outbreak occured in Germany infecting more than 3,000 people and causing more than 40 deaths. Public health officials set in motion efforts aimed at identifying the strain of *E. coli* involved and its pathogenic potential. Part of these efforts involved sequencing of DNA samples of the isolated strain on the Illumina HiSeq platform.


## Scope
In this workshop we will tackle, hands-on, the basic principles employed to generate consensus genome sequences of *E. coli* that caused the outbreak in Germany, and identify the serotypes involved in outbreaks, virulence factors, and possible anti-microbial resistance genes (AMRs). We will use the genome sequencing data of an isolate from a 16-year-old girl admitted to hospital with bloody diarrhea and abdominal pain.

## Prerequisite

This module will come after the introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment.

>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

### Set-Up
We will use the computer lab at ILRI, which is already equipped with Linux-operating desktop computers. Since we will be working from the remote servers, we will not need special setup for personal laptops. However, toward the end of the program, we can look into access to a Linux server from a Windows PC; or how to install a Linux (sub)system for any interested persons.  

### Preparations

#### ***Log into the HPC***
From the terminal (or equvalent tool) of your local computer, you can log into the HPC using the following command line, followed by pressing <ENTER>. You will be promted to type-in your password (the password will not be visible as you type it; just have faith). On a Linux system, you can use the `Ctrl-Alt-T` keyboard shortcut to open a terminal.  

```
ssh <user_name>@hpc.ilri.cgiar.org
```  

The HPC head node has 4 CPUs and we need to access more CPUs/resources in other compute nodes.
You will have to move from the cluster's master node into the node where we will be working from (it is called `compute05`). Use the following command; `-w` requests (a) specific list of host(s).  

```
interactive -w compute05 -c 2
```  

`ssh` allows you to securely connect to the remote computer over internet, while `interactive` allows you to reserve resources to work interactively in a specified node within the computing cluster using the `-w` flag.

>**Note**
>When running a job interactively, the time limit is 8 hours and Default number of CPU is 1.

#### ***Project organisation***

1. We then change into the `compute05` `scratch` directory to create our project directory. Using the option`-p`, (parent) `mkdir` will create any missing intermediate directories.
    ```
    cd /var/scratch/
    mkdir -p $USER/bacteria-wgs
    cd $USER/bacteria-wgs
    ```
2. The `databases, scripts` directories will be linked to the project directory, to limit redundancy. `-s` (soft) means that we are creating a soft link.

    ```
    ln -s /var/scratch/global/bacteria-wgs/[sd]* .
    ```
3. We will create the directories `results` and `genome` to store output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory. We will exploit the bash array data structure to create all the directories at once.
    ```
    mkdir -p genome results
    mkdir -p results/{fastqc,fastp,spades,quast,busco,prokka,bwa,blast,rgi}
    ```
4. List the contents of the `data/raw_data` directory, from where we will retrieve our ```fastq``` files.
    ```
    ls -lht data/raw_data
    
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from nucleotide sequence databases, universal `Unix` command (`wget`) can be used to download data over internet.
    ```
    # data already downloaded. Skip running these lines
    # wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292862/SRR292862_1.fastq.gz
    # wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292862/SRR292862_2.fastq.gz
    ```
2.  Download *E. coli* reference genome and the genome annotation file.

    ```
    cd genome
    ```

    We will retrieve *E. coli* reference genome and the annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type '*E. coli*' on the search box and select 'Genome' database.
    3. Right click on the genome FASTA and select 'copy link'.
    4. Use ```wget``` to fetch the files as follows:    

        ```
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
        ```
    5. Retrieve the feature annotation file GFF using ```wget``` command as follows:    

        ```
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
        ```
    6. Rename the `FASTA` and `GFF` files
        ```
        mv GCF_000005845.2_ASM584v2_genomic.fna.gz ecoli_k12_substrain_genome.fasta.gz
        mv GCF_000005845.2_ASM584v2_genomic.gff.gz ecoli_k12_substrain.gff.gz
        ```
    7. Uncompress the `gz` files
        ```
        gunzip *.gz
        ```

## Analysis

### ***Loading modules***

1. Load modules using the `module load <tool-name>` command.
    ```
    module load fastqc/0.11.7
    module load fastp/0.22.0
    ```

    **Optional**  
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 fastp/0.22.0 
    ```
2. To list the loaded modules, type the command below.
    ```
    module list
    ```
3. To unload the modules, type the command below
   ```
   module purge
   ```

### ***Prepare the reference genome***


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with `\n`. This is essential for samtools' operations.

    ```
    module load samtools/1.9
    ```
    ```
    samtools faidx ecoli_k12_substrain_genome.fasta
    ```
    The above command generates the index for reference genome with the name `ecoli_k12_substrain_genome.fasta.fai`  

    View the content of the file we just created. 
    ```
    less -S ecoli_k12_substrain_genome.fasta.fai
    ``` 
  
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 ecoli_k12_substrain_genome.fasta.fai > ecoli_k12_substrain_genome.fasta.sizes
    ```
    >**<strong style="color:magenta;opacity: 0.80;">Quiz:</strong>**
    - What is the size of the *E. coli* genome? 
    
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bwa index``` command.

    ```
    module load bwa/0.7.17
    ```
    ```
    bwa index ecoli_k12_substrain_genome.fasta

    ```

### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be critical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
    ```
    cd /var/scratch/$USER/bacteria-wgs/results/fastqc/
    ```
2. Load the `fastqc` tool
   ```
   module load fastqc/0.11.7
   ```
3. Run ```fastqc```
    ```
    fastqc \
        -t 2 \
        -o . \
        /var/scratch/$USER/bacteria-wgs/data/raw_data/SRR292862_1.fastq.gz \
        /var/scratch/$USER/bacteria-wgs/data/raw_data/SRR292862_2.fastq.gz
    ```
    ***Optional***
        Run step 3. above for the other 2 samples.

We can view the fastq report here: [SRR292862_1_fastqc.html](https://hpc.ilri.cgiar.org/~kmwangi/SRR292862_1_fastqc.html)

### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory.
    ```
    cd /var/scratch/$USER/bacteria-wgs/results/fastp/
    ```

2. Load the `fastp` tool
   ```
   module load fastp/0.22.0
   ```
3. Run ```fastp```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

    ```
    fastp \
        -w 2 \
        -i /var/scratch/$USER/bacteria-wgs/data/raw_data/SRR292862_1.fastq.gz \
        -I /var/scratch/$USER/bacteria-wgs/data/raw_data/SRR292862_2.fastq.gz \
        -o SRR292862_1.trim.fastq.gz \
        -O SRR292862_2.trim.fastq.gz \
        -h SRR292862.fastp.html \
        -j SRR292862.fastp.json \
        2> SRR292862.fastp.log
    ```
We can view the results of fastp here: [SRR292862.fastp.html](https://hpc.ilri.cgiar.org/~kmwangi/SRR292862.fastp.html)
 
### ***Perform genome assembly***  

Genome assembly refers to the process of putting back together the nucleotide sequences usually using short DNA sequences to create a representation of the original chromosome from which the sequences originated. The goal of genome assembly tools is to create long contiguous pieces of sequence (contigs) from short reads. The contigs are then ordered and oriented in relation to one another to form scaffolds. Genome assembly is a computationally intensive and difficult problem. The assembly tools have parameters that need to be tweaked and have a large effect on the outcome of any assembly. The parameters are adjusted accordingly until a desirable draft genome assembly is achieved. 

1. Change into the `spades` directory
    ```
    cd /var/scratch/$USER/bacteria-wgs/results/spades/
    ```

2. Load `Spades`
   ```
   module load spades/3.15
   ```
3. Run `Spades` to perform genome assembly
    ``` 
    spades.py -k 27 \
	-1 /var/scratch/$USER/bacteria-wgs/results/fastp/SRR292862_1.trim.fastq.gz \
	-2 /var/scratch/$USER/bacteria-wgs/results/fastp/SRR292862_2.trim.fastq.gz \
	-o ./ \
	-t 2 \
	-m 50
    ```

    View the `contigs.fasta` file

    ```
    less -S contigs.fasta
    ```
    View the top 10 headers
    ```
    grep '>' contigs.fasta | head
    ```

>**<strong style="color:magenta;opacity: 0.80;">Quiz:</strong>** 
- How many sequences do you have in your `contigs.fasta` file?
  

### ***Genome assessment***  
Genome assessment entails producing the quality metrics that gauge both the completeness and contiguity of the assembled genome. Good quality metrics ensure that we are confident in the biological insights we obtain from the assembled genome. 

1. #### ***Genome contiguity***  
    Genome contiguity is the length cutoff for the longest contigs that contain 50% of the total genome length. It is often measured as contig N50. We will use QAUST to keep track of the length of each contig and provides basic statistics. However, QUAST has many functionalities.

    ```
    module load quast/5.0.2
    ```

    ```
    quast.py \
    /var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
    -t 2 \
    -o /var/scratch/$USER/bacteria-wgs/results/quast
    ```
- Click on the following link to view the [quast report](https://hpc.ilri.cgiar.org/~kmwangi/quast/report.html).
- From the same report, we can view the statistics on a [contigs viewer using Icarus](https://hpc.ilri.cgiar.org/~kmwangi/quast/icarus_viewers/contig_size_viewer.html)

1. #### ***Genome completeness***  
    Genome completeness assesses the presence or absence of highly conserved genes (orthologs) in an assembly. This assessment is usually performed using BUSCO (Benchmarking Universal Single-Copy Orthologs). BUSCO makes use of the OrthoDB set of single-copy orthologous that are found in at least 90% of all the organisms in question. Ideally, the sequenced genome should contain most of these highly conserved genes. If your genome doesn't contain a large portion of these single-copy orthologs it may indicate that your genome is not complete. Here is BUSCO's [user guide](https://vcru.wisc.edu/simonlab/bioinformatics/programs/busco/BUSCO_v3_userguide.pdf).

    Load `BUSCO`
    ```
    module load BUSCO/5.2.2
    ```
    ```
    cd /var/scratch/$USER/bacteria-wgs/results/busco 
    ```

    ```
    busco \
    -i /var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
    -m genome \
    -o SRR292862_busco \
    -l bacteria \
    -c 2 \
    -f
    ```

    View the short summary
    ```
    less -S SRR292862_busco/short_summary.specific.bacteria_odb10.SRR292862_busco.txt

    ```

    View the full table
    ```
    less -S SRR292862_busco/run_bacteria_odb10/full_table.tsv
    ```

    List and view a amino acid of protein sequence
    ```
    ls SRR292862_busco/run_bacteria_odb10/busco_sequences/
    ```

### ***Genome annotation***  
In genome annotation, the goal is to identify and label the features of on a genome sequence.  

First, unload the modules to avoid conflict between the loaded dependencies
```
module purge
```

```
module load prokka/1.11
```
```
cd /var/scratch/$USER/bacteria-wgs/results/prokka
```  

```
prokka \
/var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
--outdir /var/scratch/$USER/bacteria-wgs/results/prokka \
--cpus 2 \
--mincontiglen 200 \
--centre C \
--locustag L \
--compliant \
--force

```

Protein abundance
```
grep -o "product=.*" L_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
```

View protein abundances
```
less -S protein_abundances.txt
```


### ***Organism identification***  
Here, we will attempt to identify the organism. Although, we know that we are working with *E. coli*, we anticipate to identify the particular _strain_ that we are working with since that is unknown at this point. We will use blastn to search a nucleotide query against a nucleotide reference database available locally.
```
module load blast/2.12.0+
```

```
cd /var/scratch/$USER/bacteria-wgs/results/blast
```

```
blastn \
-task megablast \
-query /var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
-db /export/data/bio/ncbi/blast/db/v5/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 2 \
-evalue 1e-25 \
-out /var/scratch/$USER/bacteria-wgs/results/blast/contigs.fasta.vs.nt.cul5.1e25.megablast.out

```
View the blast results
```
less -S /var/scratch/$USER/bacteria-wgs/results/blast/contigs.fasta.vs.nt.cul5.1e25.megablast.out
```

>**<strong style="color:magenta;opacity: 0.80;">Quiz:</strong>** 
- Which strain of *E. coli* was identified?


### ***Identification of virulence factors***
Virulence factors are properties that enable a microbe to establish itself within a host and contribute to its potential to cause disease (pathogenecity). We will use  the virulence factor database (VFDB) which is an integrated and comprehensive resource for curating information about virulence factors of bacterial pathogens. Briefly, we will blast our contigs against the VFDB to identify sequences carrying known virulence factors.
```
    blastn \
    -query /var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
    -db /var/scratch/${USER}/bacteria-wgs/databases/VFDB/vfdb_seta_nt \
    -out /var/scratch/$USER/bacteria-wgs/results/blast/virulence_factors_blast_VFDB.out \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -num_threads 2

```

View the blast output results
```
less -S /var/scratch/$USER/bacteria-wgs/results/blast/virulence_factors_blast_VFDB.out
```

Learn more about [blast output format 6](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

>**<strong style="color:magenta;opacity: 0.80;">Quiz:</strong>** 
Search the blast output to see if we find the following genes:
- stx1
- stx2

### ***AMR identification***  
Here, we will use [Resistance Gene Identifier (RGI)](https://card.mcmaster.ca/analyze/rgi) which uses the [Comprehensive Antibiotic Resistance Database (CARD)](https://card.mcmaster.ca/) as a reference to predict antibiotic resistome(s) from protein or nucleotide data based on homology and SNP models. 

```
module purge

module load rgi/6.0.2
```  

```
cd /var/scratch/$USER/bacteria-wgs/results/rgi

ln -s /var/scratch/global/bacteria-wgs/databases/localDB .
```


```
# Perform RGI analysis

rgi main --input_sequence /var/scratch/$USER/bacteria-wgs/results/spades/contigs.fasta \
--output_file /var/scratch/$USER/bacteria-wgs/results/rgi/SRR292862_rgi \
--local \
-a BLAST \
-g PRODIGAL \
--clean \
--low_quality \
--num_threads 1 \
--split_prodigal_jobs
	
	
# samples and AMR genes organized alphabetically:
rgi heatmap --input /var/scratch/$USER/bacteria-wgs/results/rgi \
--output /var/scratch/$USER/bacteria-wgs/results/rgi/SRR292862_rgi_alphabetic.png

```
