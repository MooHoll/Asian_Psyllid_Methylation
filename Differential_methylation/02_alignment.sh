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

if [ ! -d "lambda" ]; then                                                                                 
    mkdir lambda
fi

if [ ! -d "Dcitri" ]; then                                                                                 
    mkdir Dcitri
fi

rsync /data/ross/misc/analyses/asian_psyllid/genome/lambda.fa ./lambda
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v3.0.ref.fa ./Dcitri

echo "prep genomes"
bismark_genome_preparation ./lambda
bismark_genome_preparation ./Dcitri

echo "align to lambda read 1"
for file in $(ls *1.clean.fq.gz)
do
	base=$(basename ${file} "1.clean.fq.gz")
	bismark --multicore 8 -o ./alignment_lambda \
	--genome $SCRATCH/lambda \
	--single_end ${base}1.clean.fq.gz
done

echo "align to lambda read 2"
for file in $(ls *2.clean.fq.gz)
do
	base=$(basename ${file} "2.clean.fq.gz")
	bismark --pbat --multicore 8 -o ./alignment_lambda \
	--genome $SCRATCH/lambda \
	--single_end ${base}2.clean.fq.gz
done

echo "moving lambda outputs"
mv ./alignment_lambda /data/ross/misc/analyses/asian_psyllid/4_methylation


echo "align to Dcitri read 1"
for file in $(ls *1.clean.fq.gz)
do
	base=$(basename ${file} "1.clean.fq.gz")
	bismark --score_min L,0,-0.6 --multicore 8 -o ./alignment_Dcitri_trial \
	--genome $SCRATCH/Dcitri \
	--single_end ${base}1.clean.fq.gz
done

echo "align to Dcitri read 2"
for file in $(ls *2.clean.fq.gz)
do
	base=$(basename ${file} "2.clean.fq.gz")
	bismark --score_min L,0,-0.6 --pbat --multicore 8 -o ./alignment_Dcitri_trial \
	--genome $SCRATCH/Dcitri \
	--single_end ${base}2.clean.fq.gz
done

echo "moving Dcitri outputs"
mv ./alignment_Dcitri_trial /data/ross/misc/analyses/asian_psyllid/4_methylation

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
