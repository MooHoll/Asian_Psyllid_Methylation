#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/calc_depth_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/*.bam ./

echo "calculating coverage by scaffold"
# NOTE: taken from Kamil's pipeline: https://github.com/RossLab/PGE/blob/master/sex_chromosome_p_citry.md
for file in $(ls *sorted.bam)
do 
	base=$(basename ${file} "sorted.bam")
	samtools depth ${file} | \
    awk '/BEGIN/{scf='DC3.0sc00'; coverage_sum = 0; }{ if( scf != $1 ){ print scf "\t" coverage_sum; scf = $1; coverage_sum = $3 } else { scf = $1; coverage_sum += $3} }' > ${base}depth.txt
done

echo "moving outputs"
mv *.txt /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"