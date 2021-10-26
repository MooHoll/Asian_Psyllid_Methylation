#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/making_cov_files_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/methylation_extraction/*.cov.gz ./
rsync /data/ross/misc/analyses/asian_psyllid/logs/making_final_coverage_files.R ./

echo "doing the shiz"
gunzip *cov.gz

for file in $(ls *cov)
do
    base=$(basename ${file} "_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
    cut -f1,2,5,6 ${file} > ${base}_coverage.txt
done

R --save -q -f  making_final_coverage_files.R

echo "moving outputs"
mv ./*final* /data/ross/misc/analyses/asian_psyllid/4_methylation/methylation_extraction

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"