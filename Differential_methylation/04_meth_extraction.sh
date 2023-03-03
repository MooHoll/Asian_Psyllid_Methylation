#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/meth_extraction_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/4_methylation/alignment/*deduplicated.bam ./
rsync -r /data/ross/misc/analyses/asian_psyllid/genome/bismark_genome ./

echo "doing the shiz"
mkdir methylation_extraction

for file in $(ls *.bam)
do
    bismark_methylation_extractor -s \
    --comprehensive \
    --multicore 10 \
    --bedgraph \
    --cytosine_report \
    --genome_folder $SCRATCH/bismark_genome \
    --scaffolds \
    --output ./methylation_extraction \
    ${file}
done


echo "moving outputs"
mv ./methylation_extraction /data/ross/misc/analyses/asian_psyllid/4_methylation

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
