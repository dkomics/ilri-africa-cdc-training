---
title: README.md
tags: [ "Pathogen Genomics", "Bioinformatics", "Metadata", "Linux", "Analysis", "Activity"]
---
# **Building capacity in pathogen genomics in Africa**
---
###### ***Trainers***: [John Juma](https://github.com/ajodeh-juma), [Kennedy Mwangi](https://github.com/wanjauk) & [Gilbert Kibet](https://github.com/kibet-gilbert)
---

- [Introduction](#introduction)
- [Scope](#scope)
- [Background](#background)
- [Prerequisite](#prerequisite)
- [Set-Up](#setup)
- [Preparations](#preparations)
    - [Log into the HPC](#log-into-the-HPC)
    - [Project Organisation](#project-organisation)
    - [Data retrieval and integrity checks](#data-retrieval-and-integrity-checks)
- [Analysis](#analysis)
    - [Loading modules](#loading-modules)
    - [Prepare the reference genome](#prepare-the-reference-genome)
    - [Quality assessment](#quality-assessment)
    - [Quality and Adapter filtering](#quality-and-adapter-filtering)
    - [Consensus genome assembly](#consensus-genome-assembly)
    - [Genome assessment](#genome-assessment)
        - [Genome contiguity](#genome-contiguity)
        - [Genome completeness](#genome-completeness)
        - [Genome contamination](#genome-contamination)
    - [Genome annotation](#genome-annotation)
    - [Organism identification](#organism-identification)
        - [BLAST](#blast)
    - [Genome clean-up](#genome-clean-up)
        - [Read mapping](#read-mapping)
        - [Construct a coverage table](#Construct-a-coverage-table)
        - [Non-target contig removal](#Non-target-contig-removal)
        - [Filter the genome assembly by length](#Filter-the-genome-assembly-by-length)
        - [Filter the genome assembly by coverage](#Filter-the-genome-assembly-by-coverage)
        - [Construct a list of contigs to keep](#Construct-a-list-of-contigs-to-keep)
        - [Filter your assembly based on a list of contigs](#Filter-your-assembly-based-on-a-list-of-contigs)
        - [Check for contamination using UniVec](#Check-for-contamination-using-UniVec)



## Introduction


## Scope
In this workshop we will tackle, hands-on, the basic principles employed by the numerous bioinformatic pipelines: to generate consensus genome sequences of *E. coli* and identify variants using an actual dataset generated in our facility.

> **Note**

> This is part of the initiative fronted by the [Africa CDC](https://africacdc.org/) with generous support from the [Rockeffeler foundation](https://www.rockefellerfoundation.org/) to build capacity in pathogen genomics in Africa.


## Background



## Prerequisite

This module will come after the introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment.

>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned, for example `Bio4Info$$`. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

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
    mkdir -p $USER/AfricaCDC_training
    cd $USER/AfricaCDC_training
    ```
2. The `<directories to soft-link>` directories will be linked to the project directory, to limit redundancy. `-s` (soft) means that we are creating a soft link.
    ```
    ln -s /var/scratch/global/AfricaCDC_training/[adps]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw data in ```fastq``` format, output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory. We will exploit the bash array data structure to create all the directories at once.
    ```
    mkdir data genome results
    mkdir -p results/{fastqc,fastp,kraken,samtools,ivar,snpeff,pangolin,nextclade,multiqc,bowtie2,bedtools}
    ```
4. Change into the `data` directory, from where we will retrieve our ```fastq``` files.
    ```
    cd data
    ls
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from nucleotide sequence databases, universal `Unix` command (`wget`) can be used to download data over internet.
    ```
    wget --no-check-certificate <path to compressed data>
    ```
2. After downloading your data, say from a sequencing facility site, it is often good practice to verify that your data was not intentionally/accidentally tampered with. To do this, your data service provider will likely accompany your data with a file containing a verification code: `checksum_file` (***will be provided***). The `md5sum` command, using the `-c` (check) tag, allows for checking the integrity of a file downloaded or acquired from a different source.
    ```
    wget --no-check-certificate <path to md5 checksum file of the compressed data>
    ls
    md5sum -c <md5 checksum file of the compressed data>
    ```
3. Next, we will unzip the file using `tar` with the `-xf` (extract, file; respectively) tags, which tells `tar` extract the given file.
    ```
    tar -xf <compressed data>
    ls
    ```
4.  Download *E. coli* reference genome and the genome annotation file.

    We will retrieve *E. coli* reference genome and the annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type '*E. coli*' on the search box and select 'Genome' database.
    3. Select the Genbank hyperlink.
    4. Select the genome version.
    5. Right click on the genome FASTA and select 'copy link'.
    6. Change into the ```genome``` directory using the command
    ```cd ../genome```.
    7. Use ```wget``` to fetch the file.
    8. Retrieve the feature annotation file GFF using ```wget``` command.
    9. Dowload the md5checksum using `wget` command and check for integrity of your reference genome (FASTA) and annotation (GFF) files.

        ```
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.fna.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.fna.gz" | md5sum -c -
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.gff.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.gff.gz" | md5sum -c -
        ```

    11. If integrity check of the files has passed (`OK`), Uncompress the ```.gz``` files
       ```
       gunzip *.gz
       ```
    12. Rename the `FASTA` and `GFF` files
        ```
        mv 
        mv 
        ```

## Analysis

#### ***Loading modules***
1. Clear the environment.
    ```
    module purge
    ```
2. Load modules using the `module load <tool-name>`command.
    ```
    module load fastqc/0.11.7
    module load fastp/0.22.0
    module load kraken/2.0.8-beta
    module load bowtie2/2.3.4.1
    module load samtools/1.11
    module load ivar/1.3.1
    module load bedtools/2.29.0
    module load R/3.6
    module load bcftools/1.11
    module load snpeff/4.1g
    module load multiqc/1.12
    module load nextclade/1.11.0
    module load python/3.9
    ```

    **Optional**
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 fastp/0.22.0 \
    kraken/2.0.8-beta bowtie2/2.3.4.1 samtools/1.11 ivar/1.3.1 \
    bedtools/2.29.0 R/3.6 bcftools/1.11 snpeff/4.1g multiqc/1.12 \
    nextclade/1.11.0 python/3.9
    ```
3. To list the loaded modules, type the below command.
    ```
    module list
    ```

#### ***Prepare the reference genome***


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for samtools' operations.

    ```
    samtools faidx <E-coli>.fasta
    ```
    The above command generates the index for reference genome with the name `<E-coli>.fasta.fai`.
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 <E-coli>.fasta.fai > <E-coli>.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bowtie2-build``` command.

    ```
    mkdir /var/scratch/$USER/AfricaCDC_training/genome/bowtie2
    ```

    ```
    bowtie2-build \
      --threads 1 \
      /var/scratch/$USER/AfricaCDC_training/genome/<E-coli>.fasta \
      /var/scratch/$USER/AfricaCDC_training/genome/bowtie2/<E-coli>
    ```
    The above command generates index files with the suffix `.bt2` for the reference genome with the prefix `<E-coli>.`
4. Build SnpEff database for the reference genome

    [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/), a variant annotation and predictor needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes, so chances are that your organism of choice already has a SnpEff database available.

    >**Note** We will use pre-built *E. coli* SnpEff database


    ***Optional***
    In the (unlikely?) event that you need to build one yourself, you can build one using the commands found [here](http://pcingola.github.io/SnpEff/se_buildingdb/)

#### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastqc/
    ```
2. Run ```fastqc```
    ```
    fastqc \
        -t 1 \
        -o . \
        /var/scratch/$USER/AfricaCDC_training/data/E-coli_R1.fastq.gz \
        /var/scratch/$USER/AfricaCDC_training/data/E-coli_R2.fastq.gz
    ```
    ***Optional***
        Run step 3. above for the other 2 samples.

#### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory.
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastp/
    ```

2. Run ```fastp```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

    ```
    fastp \
        -w 1 \
        -i /var/scratch/$USER/AfricaCDC_training/data/E-coli_R1.fastq.gz \
        -I /var/scratch/$USER/AfricaCDC_training/data/E-coli_R2.fastq.gz \
        -o E-coli_R1.trim.fastq.gz \
        -O E-coli_R2.trim.fastq.gz \
        -h E-coli.fastp.html \
        -j E-coli.fastp.json \
        2> E-coli.fastp.log
    ```

    ***Optional***
        Run steps 3 and 4 above for the other 2 samples.

#### ***Consensus genome assembly***
To generate a consensus sequence iVar uses the output of samtools mpileup command. The mpileup output must be piped into ivar consensus. There are five parameters that can be set:
- minimum quality ```-q``` (Default: 20).
- minimum frequency threshold ```-t``` (Default: 0).
- minimum depth to call a consensus ```-m``` (Default: 10).
- a flag ```-n``` to exclude nucleotides from regions with depth less than the minimum depth and a character to call in regions with coverage lower than the speicifed minimum depth (Default: 'N').

Minimum quality is the minimum quality of a base to be considered in calculations of variant frequencies at a given position. Minimum frequency threshold is the minimum frequency that a base must match to be called as the consensus base at a position. If one base is not enough to match a given frequency, then an ambigious nucleotide is called at that position. Minimum depth is the minimum required depth to call a consensus. If ```-k``` flag is set then these regions are not included in the consensus sequence. If ```-k``` is not set then by default, a 'N' is called in these regions. You can also specfy which character you want to add to the consensus to cover regions with depth less than the minimum depth. This can be done using ```-n``` option. It takes one of two values: ```-``` or ```N```.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Generate pileup and consensus genome sequences

    ```
    samtools \
            mpileup \
            --reference /var/scratch/$USER/AfricaCDC_training/genome/<E-coli>.fasta \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            -aa \
            /var/scratch/$USER/AfricaCDC_training/results/ivar/E-coli.primertrimmed.sorted.bam \
            | tee E-coli.mpileup \
            | ivar \
                consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p E-coli.cons
    ```
The ```tee``` command reads from the standard input and writes to both standard output and one or more files at the same time. ```tee``` is mostly used in combination with other commands through piping.

3. Rename the consensus genome header

    ```
    sed -i '/^>/s/Consensus_\(.*\)_threshold.*/\1/' E-coli.cons.fa
    ```


#### ***Genome assessment***


##### ***Genome contiguity***


##### ***Genome completeness***


##### ***Genome contamination***


#### ***Genome annotation***


#### ***Organism identification***


##### ***BLAST***


#### ***Genome clean-up***


##### ***Read mapping***


##### ***Construct a coverage table***


##### ***Non-target contig removal***


##### ***Filter the genome assembly by length***


##### ***Filter the genome assembly by coverage***


##### ***Construct a list of contigs to keep***


##### ***Filter your assembly based on a list of contigs***


##### ***Check for contamination using UniVec***




