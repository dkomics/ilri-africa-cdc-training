#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J prokka
#SBATCH -n 8
#SBATCH -o slurm-%x.out
#SBATCH -e slurm-%x.err

# load the required module
#module load quast/5.0.2
#module load BUSCO/5.2.2
module load prokka/1.14.6


res1=$(date +%s.%N)

# Run quast
#quast.py \
#/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
#-t 8 \
#-o /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/quast


# Run BUSCO
# Note: to use the full path fr output directory, one has to edit and place it in the congig.ini file 
#busco \
#-i /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
#-m genome \
#-o busco-results \
#-l bacteria \
#-c 8

# Run PROKKA
prokka \
/var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/spades/contigs.fasta \
--outdir /var/scratch/global/${USER}/ilri-africa-cdc-training/bacteria/data/prokka \
--cpus 8\ 
--mincontiglen 200 \
--force






res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

echo
echo
echo
LC_NUMERIC=C printf "Total runtime: %d:%02d:%02d:%07.4f\n" $dd $dh $dm $ds
