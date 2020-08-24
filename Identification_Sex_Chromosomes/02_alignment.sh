#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/alignment_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
# Alignment to reference genome version 2
# Download genome from citrusgreening.org and cp to server

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v2.0.ref.fasta ./
rsync /data/ross/sequencing/raw/asian_psyllid/*clean.fq.gz ./

echo "making the genome index"
bowtie2-build ./Diaci_v2.0.ref.fasta d_citri

echo "starting alignment"
for file in $(ls *_1.clean.fq.gz)
do
	base=$(basename $file "_1.clean.fq.gz")
    bowtie2 --sensitive --threads 8 \
    -x d_citri \
    -1 ${base}_1.clean.fq.gz -2 ${base}_2.clean.fq.gz \
    -S ${base}.sam
done

echo "convert sams to bams"
for file in $(ls *.sam)
do
	base=$(basename $file ".sam")
    samtools view -bS ${base}.sam > ${base}.bam
done

echo "moving outputs"
mv *.bam /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"