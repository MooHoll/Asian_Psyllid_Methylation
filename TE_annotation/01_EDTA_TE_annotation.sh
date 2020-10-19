#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/te_annotation_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`
#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v3.0.ref.fa ./
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitr_OGSv3.0_beta_cds.fa ./
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitri_genes.bed ./

echo "running EDTA"

perl /data/ross/misc/analyses/asian_psyllid/te_analysis/EDTA/EDTA_raw.pl \
--genome Diaci_v3.0.ref.fa \
--cds Dcitr_OGSv3.0_beta_cds.fa \
--anno 1 \
--exclude Dcitri_genes.bed \
--threads 32

echo "moving outputs"
mv * /data/ross/misc/analyses/asian_psyllid/te_analysis

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"