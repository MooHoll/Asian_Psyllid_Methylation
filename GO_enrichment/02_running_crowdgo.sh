#!/bin/bash

#PBS -N GO_annotation
#PBS -l walltime=00:05:00
#PBS -l vmem=100gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=24
#PBS -q devel

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load python/gcc/3.5.3
module load interproscan/5.46.81
module load snakemake/4.3.1
module load hmmer/3.1b2
source /scratch/monoallelic/hm257/daphnia/bin/crowdgo_env/bin/activate

snakemake -s CrowdGO.snakefile