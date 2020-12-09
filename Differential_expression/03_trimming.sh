#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/trimming_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER/$JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/raw_data/RNA-Seq/*.fq.gz ./

echo "doing the shiz"
for file in $(ls *1.clean.fq.gz)
do
	base=$(basename $file "_1.clean.fq.gz")
	cutadapt -j 10 \
	-u 10 -U 10 \
	-o trim_${base}1.fq.gz \
	-p trim_${base}2.fq.gz \
	${base}1.fq.gz \
	${base}2.fq.gz
done

echo "moving outputs"
mv trim_* /data/ross/misc/analyses/asian_psyllid/raw_data/trimmed_RNA-Seq

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
