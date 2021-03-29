#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/cpg_annotation_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------
# conda activate R_env

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/total_cpgs_in_genome.txt ./
rsync /data/ross/misc/analyses/asian_psyllid/windows.bed ./
rsync /data/ross/misc/analyses/asian_psyllid/logs/total_cpgs_per_annotation.R ./

# Split windows file into 10 files making sure not to break lines
split -n l/20 windows.bed

echo "doing the shiz"
R --save -f total_cpgs_per_annotation.R


echo "moving outputs"
mv windows_with_total_cpgs.txt /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"
