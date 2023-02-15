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
    - [Genome visualization](#genome-visualization)
    - [AMR identification](#AMR-identification)


## Introduction


## Scope
In this workshop we will tackle, hands-on, the basic principles employed by the numerous bioinformatic pipelines: to generate consensus genome sequences of *E. coli* and identify variants using an actual dataset generated in our facility.

> **Note**

> This is part of the initiative fronted by the [Africa CDC](https://africacdc.org/) with generous support from the <add org here> to build capacity in pathogen genomics in Africa.


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
    mkdir -p $USER/ilri-africa-cdc-training
    cd $USER/ilri-africa-cdc-training
    ```
2. The `<directories to soft-link>` directories will be linked to the project directory, to limit redundancy. `-s` (soft) means that we are creating a soft link.

    ```
    ln -s /var/scratch/global/ilri-africa-cdc-training/[sd]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw data in ```fastq``` format, output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory. We will exploit the bash array data structure to create all the directories at once.
    ```
    mkdir -p data genome results
    mkdir -p results/{fastqc,fastp,spades,quast,busco,prokka,bwa,rgi}
    ```
4. Change into the `data` directory, from where we will retrieve our ```fastq``` files.
    ```
    cd data/raw_data
    ls
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from nucleotide sequence databases, universal `Unix` command (`wget`) can be used to download data over internet.
    ```
    # data already downloaded. Skip running these lines
    #wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_1.fastq.gz
    #wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR292/SRR292770/SRR292770_2.fastq.gz
    ```
2.  Download *E. coli* reference genome and the genome annotation file.

    We will retrieve *E. coli* reference genome and the annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type '*E. coli*' on the search box and select 'Genome' database.
    3. Right click on the genome FASTA and select 'copy link'.
    4. Change into the ```genome``` directory using the command
    ```cd ../../genome```.
    5. Use ```wget``` to fetch the files as follows:    

            ```
            wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
           
            ```
    6. Retrieve the feature annotation file GFF using ```wget``` command as follows:    

            ```
             wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
            ```
    7. Rename the `FASTA` and `GFF` files
        ```
        mv GCF_000005845.2_ASM584v2_genomic.fna.gz E-coli-genome.fasta.gz
        mv GCF_000005845.2_ASM584v2_genomic.gff.gz E-coli.gff.gz
        ```
    8. Uncompress the `gz` files
        ```
        gunzip *.gz
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
    
    ```

    **Optional**
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 fastp/0.22.0 
    ```
3. To list the loaded modules, type the below command.
    ```
    module list
    ```

#### ***Prepare the reference genome***


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for samtools' operations.

    ```
    module load samtools/1.9

    samtools faidx E-coli-genome.fasta
    ```
    The above command generates the index for reference genome with the name `E-coli-genome.fasta.fai`.
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 E-coli-genome.fasta.fai > E-coli-genome.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bwa index``` command.

    ```
    mkdir /var/scratch/$USER/ilri-africa-cdc-training/genome/bwa
    ```

    ```
    bwa index E-coli-genome.fasta

    ```

#### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
    ```
    cd /var/scratch/$USER/ilri-africa-cdc-training/results/fastqc/
    ```
2. Run ```fastqc```
    ```
    fastqc \
        -t 1 \
        -o . \
        /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_1.fastq.gz \
        /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_2.fastq.gz
    ```
    ***Optional***
        Run step 3. above for the other 2 samples.

#### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory.
    ```
    cd /var/scratch/$USER/ilri-africa-cdc-training/results/fastp/
    ```

2. Run ```fastp```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

    ```
    fastp \
        -w 1 \
        -i /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_1.fastq.gz \
        -I /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_2.fastq.gz \
        -o SRR292862_1.trim.fastq.gz \
        -O SRR292862_2.trim.fastq.gz \
        -h SRR292862.fastp.html \
        -j SRR292862.fastp.json \
        2> SRR292862.fastp.log
    ```

    ***Optional***
        Run steps 3 and 4 above for the other 2 samples.

#### ***Consensus genome assembly***  

Generate consensus genome sequences  

    ```
    module load spades/3.15
    ```

1. Change into the `spades` directory
    ```
    cd /var/scratch/$USER/ilri-africa-cdc-training/results/spades/
    ```

2. Run `Spades` to perform genome assembly
    ``` 
    spades.py -k 27 \
	-1 /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_1.trim.fastq.gz \
	-2 /var/scratch/$USER/ilri-africa-cdc-training/data/raw_data/SRR292862_2.trim.fastq.gz \
	-o ./ \
	-t 8 \
	-m 384
    ```



#### ***Genome assessment***  

1. ##### ***Genome contiguity***  

    ```
    module load quast/5.0.2
    ```

    ```
    quast.py \
    /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta \
    -t 8 \
    -o /var/scratch/${USER}/ilri-africa-cdc-training/results/quast
    ```

2. ##### ***Genome completeness***  

    ```
    module load BUSCO/5.2.2
    ```

    ```
    busco \
    -i /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta \
    -m genome \
    -o /var/scratch/${USER}/ilri-africa-cdc-training/results/busco \
    -l bacteria \
    -c 8
    ```

3. ##### ***Genome contamination***  


#### ***Genome annotation***  

    ```
    module load prokka/1.11
    ```

    ```
    prokka \
    /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta \
    --outdir /var/scratch/${USER}/ilri-africa-cdc-training/results/prokka \
    --cpus 8 \
    --mincontiglen 200 \
    --centre XXX \
    --force
    ```

#### ***Organism identification***  

    ```
    module load blast/2.12.0+
    ```

    ```
    bash /var/scratch/kmwangi/ilri-africa-cdc-training/scripts/blob_blast.sh \
    /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta
    ```

1. ##### ***BLAST***  

    ```
    # run the script, note that it will automatically use nohup since it will take about 30 minutes to run
    bash /var/scratch/kmwangi/ilri-africa-cdc-training/scripts/blob_blast.sh \
    /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta

    # view the reuslts, the last column is the species identification
    tabview /var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta.vs.nt.cul5.1e5.megablast.out
    ```

#### ***Genome clean-up***  

    ```
    module load bwa/0.7.4
    module load samtools/1.9
    ```

1. ##### ***Read mapping***  
    ```
    fasta=/var/scratch/kmwangi/ilri-africa-cdc-training/genome/E-coli-genome.fasta
    forward=/var/scratch/gkibet/ilri-africa-cdc-training/results/fastp/SRR292862_1.trim.fastq.gz
    reverse=/var/scratch/gkibet/ilri-africa-cdc-training/results/fastp/SRR292862_1.trim.fastq.gz

    # Step 1: Index your reference genome. This is a requirement before read mapping.
    bwa index $fasta

    # Step 2: Map the reads and construct a SAM file.
    bwa mem -t 8 $fasta $forward $reverse > /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/raw_mapped.sam

    # Step 3: Remove sequencing reads that did not match to the assembly and convert the SAM to a BAM.
    samtools view -@ 8 -Sb  /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/raw_mapped.sam \
    | samtools sort -@ 8 -o /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/sorted_mapped.bam

    # Step 4: Index the new bam file
    samtools index /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/sorted_mapped.bam
    ```

2. ##### ***Construct a coverage table***  

    ```
    module load bedtools/2.29.0
    ```

    ```
    bedtools genomecov \
    -ibam /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/sorted_mapped.bam \
    > /var/scratch/${USER}/ilri-africa-cdc-training/results/bwa/coverage.out
    ```

3. ##### ***Non-target contig removal***  


4. ##### ***Filter the genome assembly by length***  


5. ##### ***Filter the genome assembly by coverage***  


6. ##### ***Construct a list of contigs to keep***  


7. ##### ***Filter your assembly based on a list of contigs***  


8. ##### ***Check for contamination using UniVec***  


#### ***Genome visualization***  

#### ***AMR identification***  

```
module purge

module load rgi/6.0.2

```  

```
    # Perform RGI analysis

    rgi_output=/var/scratch/${USER}/ilri-africa-cdc-training/results/rgi
    contigs_file=/var/scratch/${USER}/ilri-africa-cdc-training/results/spades/contigs.fasta

    rgi main --input_sequence $contigs_file \
        --output_file $rgi_output \
	    --local \
	    -a BLAST \
	    -g PRODIGAL \
	    --clean \
	    --low_quality \
	    --num_threads 6 \
	    --split_prodigal_jobs
	
	
    # samples and AMR genes organized alphabetically:
    rgi heatmap --input $rgi_output \
        --output $rgi_output/rgi_alphabetic.png

```
