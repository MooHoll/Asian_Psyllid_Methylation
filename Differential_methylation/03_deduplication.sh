#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/deduplication_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment/*.bam ./

echo "deduplicating now"
for file in $(ls *.bam)
do
	deduplicate_bismark --bam -s ${file}
done

echo "moving outputs"
mv ./*dedup* /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
