#!/bin/bash

#PBS -N blast
#PBS -l walltime=00:05:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed
module load blast+/2.9.0

blastp \
-query Dcitr_OGSv3.0_beta_pep.fa \
-db "swissprot" \
-remote \
-max_target_seqs 3 \
-outfmt 5 \
-out d_citri.xml \
-evalue 1e-3

# also check for dsx and fru in the D. citri genome
# made fasta file of all Drosophila melanogastor dsx, fru and tra isoforms from NCBI

makeblastdb -in dsx_drosophila.fa -parse_seqids -dbtype prot
blastp -query Dcitr_OGSv3.0_beta_pep.fa \
-db "dsx_drosophila.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_dsx.txt
# 11 total
# Dcitr02g06160.1.1 x3 
# Dcitr02g06180.1.1 x1 
# Dcitr02g07270.1.1 x6
# Dcitr06g02910.1.1 x1 (low score)

makeblastdb -in fru_drosophila.fa -parse_seqids -dbtype prot
blastp -query Dcitr_OGSv3.0_beta_pep.fa \
-db "fru_drosophila.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_fru.txt
# Loads (211)

makeblastdb -in tra_drosophila.fa -parse_seqids -dbtype prot
blastp -query Dcitr_OGSv3.0_beta_pep.fa \
-db "tra_drosophila.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_tra.txt
# 0

# Add in the reciprocal

makeblastdb -in Dcitr_OGSv3.0_beta_pep.fa -parse_seqids -dbtype prot
blastp -query dsx_drosophila.fa \
-db "Dcitr_OGSv3.0_beta_pep.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_reciprocal_dsx.txt
# 6 total
# Dcitr03g16970.1.1 x6 (NOTE: this gene has two annotated isoforms)

blastp -query fru_drosophila.fa \
-db "Dcitr_OGSv3.0_beta_pep.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_reciprocal_fru.txt
# 28 total
# Dcitr01g04580.1.1 x 28 (NOTE: this gene has one annotated isoform)

blastp -query tra_drosophila.fa \
-db "Dcitr_OGSv3.0_beta_pep.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out d_citri_reciprocal_tra.txt
# 0