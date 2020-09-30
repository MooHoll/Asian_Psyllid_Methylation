#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/blast_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitr_OGSv3.0_beta_singlecopy.pep.fa ./
rsync /data/ross/misc/analyses/asian_psyllid/sex_chromosome_analysis/x_chrom_synteny/Pachypsylla_venusta_singlecopy.pep.fa ./

echo "making blast database"
makeblastdb -in Pachypsylla_venusta_singlecopy.pep.fa -dbtype prot -out Pachypsylla_venusta_protein_blastdb

echo "running blast"
# following parameters set in: https://doi.org/10.1101/2020.03.24.006411
blastp -db Pachypsylla_venusta_protein_blastdb \
-query Dcitr_OGSv3.0_beta_singlecopy.pep.fa  \
-evalue 1e-10 \
-outfmt 6 \
-num_alignments 5 \
-num_threads 10 \
-out Dcitri_to_Pvenusta.blast

echo "moving outputs"
mv * /data/ross/misc/analyses/asian_psyllid/sex_chromosome_analysis/x_chrom_synteny

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"