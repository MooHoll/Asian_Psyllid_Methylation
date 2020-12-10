#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/alignment_WGBS_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/raw_data/WGBS/*.fq.gz ./

mkdir lambda
mkdir Dcitri

rsync /data/ross/misc/analyses/asian_psyllid/genome/lambda.fa ./lambda
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v3.0.ref.fa ./Dcitri

echo "prep genomes"
bismark_genome_preparation ./lambda
bismark_genome_preparation ./Dcitri

echo "align to lambda"
for file in $(ls *1.clean.fq.gz)
do
	base=$(basename ${file} "1.clean.fq.gz")
	bismark --multicore 8 -o ./alignment_lambda \
	--genome $SCRATCH/lambda \
	-1 ${base}1.clean.fq.gz \
	-2 ${base}1.clean.fq.gz
done

echo "moving lambda outputs"
mv ./alignment_lambda /data/ross/misc/analyses/asian_psyllid/4_methylation

echo "align to Dcitri"
for file in $(ls *1.clean.fq.gz)
do
	base=$(basename ${file} "1.clean.fq.gz")
	bismark --multicore 8 -o ./alignment_Dcitri \
	--genome $SCRATCH/Dcitri \
	-1 ${base}1.clean.fq.gz \
	-2 ${base}1.clean.fq.gz
done

echo "moving Dcitri outputs"
mv ./alignment_Dcitri /data/ross/misc/analyses/asian_psyllid/4_methylation

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
