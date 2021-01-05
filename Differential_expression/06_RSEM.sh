#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/rsem_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v3.0.ref.fa ./
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitr_OGSv3.0_beta.gff3 ./
rsync /data/ross/misc/analyses/asian_psyllid/raw_data/trimmed_RNA-Seq/*.fq.gz ./
#---------------------------------------------
# https://github.com/deweylab/RSEM

echo "prepare rsem reference"
rsem-prepare-reference --gff3 Dcitr_OGSv3.0_beta.gff3 --star -p 32 Diaci_v3.0.ref.fa Dcitri

echo "run RSEM"
for file in $(ls *1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    rsem-calculate-expression --star -p 20 --paired-end ${base}_1.fq.gz ${base}_2.fq.gz Dcitri ${base}
done

#---------------------------------------------

echo "moving outputs"
mv * /data/ross/misc/analyses/asian_psyllid/3_expression/rsem

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
