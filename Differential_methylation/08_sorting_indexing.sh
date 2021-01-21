#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/sorting_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment/*deduplicated.bam ./

echo "doing the shiz"

for file in $(ls *.bam)
do
    base=$(basename $file ".clean_bismark_bt2.deduplicated.bam")
    /ceph/software/samtools/samtools-1.7/samtools sort -o ${base}_deduplicated_sorted.bam -@ 10 ${file}
done

wait

for file in $(ls *sorted.bam)
do
    /ceph/software/samtools/samtools-1.7/samtools index ${file}
done

echo "moving outputs"
mv ./*sorted* /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
