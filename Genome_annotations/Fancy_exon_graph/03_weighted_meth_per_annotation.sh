#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/weighted_meth_per_window_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
# conda activate R_env

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/coverage_files/*.txt ./
rsync /data/ross/misc/analyses/asian_psyllid/windows_with_total_cpgs.txt ./
rsync /data/ross/misc/analyses/asian_psyllid/logs/weighted_meth_per_sample.R ./

echo "doing the shiz"
R --no-restore --save -f weighted_meth_per_sample.R

echo "moving outputs"
mv ./*weighted_meth.txt /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
