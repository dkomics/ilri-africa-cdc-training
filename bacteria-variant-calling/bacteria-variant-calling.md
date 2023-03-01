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
  - [Scope](#scope)
  - [Background](#background)
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



## Introduction


## Scope
In this workshop we will tackle, hands-on, the basic principles employed by the numerous bioinformatic pipelines: to generate consensus genome sequences of *E. coli* and identify variants using an actual dataset generated in our facility.

> **Note**

> This is part of the initiative fronted by the [Africa CDC](https://africacdc.org/) with generous support from the <add org here> to build capacity in pathogen genomics in Africa.


## Background



## Prerequisite

This module will come after the introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment.

>**Note**

>Once inside the `hpc`, all instances of ```${USER}``` will be equivalent to the hpc username that you were assigned, for example `Bio4Info$$`. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

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
interactive -w compute05
```  

`ssh` allows you to securely connect to the remote computer over internet, while `interactive` allows you to reserve resources to work interactively in a specified node within the computing cluster using the `-w` flag.

>**Note**
>When running a job interactively, the time limit is 8 hours and Default number of CPU is 1.

#### ***Project organisation***

1. We then change into the `compute05` `scratch` directory to create our project directory. Using the option`-p` (parent) `mkdir` will create any missing intermediate directories.
    ```
    cd /var/scratch/
    mkdir -p ${USER}/bacteria-variant-calling
    cd ${USER}/bacteria-variant-calling
    ```
2. The `<directories to soft-link>` directories will be linked to the project directory, to limit redundancy. `-s` (soft) means that we are creating a soft link.

    ```
    ln -s /var/scratch/global/bacteria-variant-calling/[d]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw data in ```fastq``` format, output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory. We will exploit the bash array data structure to create all the directories at once.
    ```
    mkdir -p genome results
    mkdir -p results/{fastqc,trimmomatic,bwa,bcf,vcf}
    ```
4. List the contents of the `data/raw_data` directory, from where we will retrieve our ```fastq``` files.
    ```
    ls data/raw_data
    
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from nucleotide sequence databases, universal `Unix` command (`wget`) can be used to download data over internet.
    ```
    # data already downloaded. Skip running these lines
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584863/SRR2584863_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584863/SRR2584863_2.fastq.gz
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
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
        ```
    5. Retrieve the feature annotation file GFF using ```wget``` command as follows:    

        ```
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.gff.gz
        ```
    6. Rename the `FASTA` and `GFF` files
        ```
        mv GCA_000017985.1_ASM1798v1_genomic.fna.gz ecoli_rel606.fasta.gz
        mv GCA_000017985.1_ASM1798v1_genomic.gff.gz ecoli_rel606.gff.gz
        ```
    7. Uncompress the `gz` files
        ```
        gunzip *.gz
        ```

## Analysis

#### ***Loading modules***

1. Load modules using the `module load <tool-name>`command.
    ```
    module load fastqc/0.11.7
    module load trimmomatic/0.39
    
    ```

    **Optional**
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 trimmomatic/0.39 
    ```
2. To list the loaded modules, type the below command.
    ```
    module list
    ```

#### ***Prepare the reference genome***


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for samtools' operations.

    ```
    module load samtools/1.9

    samtools faidx ecoli_rel606.fasta
    ```
    The above command generates the index for reference genome with the name `ecoli_rel606.fasta.fai`.
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 ecoli_rel606.fasta.fai > ecoli_rel606.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bwa index``` command.

    ```
    module load bwa/0.7.17
    
    bwa index ecoli_rel606.fasta

    ```

#### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
```
    cd /var/scratch/${USER}/bacteria-variant-calling/results/fastqc/
```
2. Run ```fastqc```
```
    fastqc \
        -t 1 \
        -o . \
        /var/scratch/${USER}/bacteria-variant-calling/data/raw_data/SRR2584863_1.fastq.gz \
        /var/scratch/${USER}/bacteria-variant-calling/data/raw_data/SRR2584863_2.fastq.gz

```
    ***Optional***
        Run step 3. above for the other 2 samples.

#### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```trimmomatic``` directory.
```
    cd /var/scratch/${USER}/bacteria-variant-calling/results/trimmomatic/
```

2. Run ```trimmomatic```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

```

    cp /export/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa .

    trimmomatic PE \
    /var/scratch/${USER}/bacteria-variant-calling/data/raw_data/SRR2584863_1.fastq.gz \
    /var/scratch/${USER}/bacteria-variant-calling/data/raw_data/SRR2584863_2.fastq.gz \
    SRR2584863_1.trim.fastq.gz \
    SRR2584863_1un.trim.fastq.gz \
    SRR2584863_2.trim.fastq.gz \
    SRR2584863_2un.trim.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15


```

 #### ***Align reads to reference genome***   

Change directory to the `bwa` results folder

```
    cd /var/scratch/${USER}/bacteria-variant-calling/results/bwa/

```

```
    bwa mem \
    /var/scratch/${USER}/bacteria-variant-calling/genome/ecoli_rel606.fasta \
    /var/scratch/${USER}/bacteria-variant-calling/results/trimmomatic/SRR2584863_1.trim.fastq.gz \
    /var/scratch/${USER}/bacteria-variant-calling/results/trimmomatic/SRR2584863_2.trim.fastq.gz \
    > SRR2584863.aligned.sam

```

Convert the SAM file to BAM format

```
samtools view -S -b \
SRR2584863.aligned.sam \
> SRR2584863.aligned.bam


```

Sort the bam file by coordinates

```
    samtools sort \
    -o SRR2584863.aligned.sorted.bam \
    SRR2584863.aligned.bam 

```


Alignment statistics

```
    samtools flagstat SRR2584863.aligned.sorted.bam

```


### Variant calling  

Change directory to the `vcf` results folder

```
    cd /var/scratch/${USER}/bacteria-variant-calling/results/vcf/

```

Step 1: Calculate the read coverage of positions in the genome  

```
    module load bcftools/1.15.1

    bcftools mpileup -O b \
    -o SRR2584863_raw.bcf \
    -f /var/scratch/${USER}/bacteria-variant-calling/genome/ecoli_rel606.fasta \
    /var/scratch/${USER}/bacteria-variant-calling/results/bwa/SRR2584863.aligned.sorted.bam 

```

Step 2: Detect the single nucleotide variants (SNVs)  

```
    bcftools call \
    --ploidy 1 \
    -m \
    -v \
    -o SRR2584863_variants.vcf \
    SRR2584863_raw.bcf 

```

Step 3: Filter and report the SNV variants in variant calling format (VCF)  

```
    vcfutils.pl varFilter \
    SRR2584863_variants.vcf \
    > SRR2584863_final_variants.vcf

```

Explore the VCF format:  

```
     less -S SRR2584863_final_variants.vcf

```

Assess the alignment (visualization)

Using `tview`.  

```  
    samtools index \
    /var/scratch/${USER}/bacteria-variant-calling/results/bwa/SRR2584863.aligned.sorted.bam

    samtools tview \
    /var/scratch/${USER}/bacteria-variant-calling/results/bwa/SRR2584863.aligned.sorted.bam \
    /var/scratch/${USER}/bacteria-variant-calling/genome/ecoli_rel606.fasta

```  

Using `igv`.  
Install IGV

