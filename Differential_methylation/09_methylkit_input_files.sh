#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/methylkit_inputs_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment/*sorted.bam ./
rsync /data/ross/misc/analyses/asian_psyllid/logs/making_files.R ./

for file in $(ls *_1_deduplicated_sorted.bam)
do
    base=$(basename ${file} "_1_deduplicated_sorted.bam")
    samtools merge -f ${base}.bam ${base}_1_deduplicated_sorted.bam ${base}_2_deduplicated_sorted.bam
done

echo "doing the shiz"
R --save -f making_files.R

echo "moving outputs"
mv ./*.txt /data/ross/misc/analyses/asian_psyllid/4_methylation/methylkit

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
