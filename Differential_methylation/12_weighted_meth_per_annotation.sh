#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/weighted_meth_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/coverage_files/*.txt ./
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitr_OGSv3.0_beta_with_total_cpgs.txt ./
rsync /data/ross/misc/analyses/asian_psyllid/logs/weighted_meth_per_sample.R ./

echo "doing the shiz"
R --save -f weighted_meth_per_sample.R

echo "moving outputs"
mv ./*weighted_meth.txt /data/ross/misc/analyses/asian_psyllid/4_methylation

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
