#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/calc_depth_for_Y$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying data in"
rsync /data/ross/misc/analyses/asian_psyllid/*.bam ./
rsync /data/ross/misc/analyses/asian_psyllid/Dcitri_chroms_length_tabfile.txt ./

# The -g input files is a tab sep txt file with chr and chr_size columns (no headers)
bedtools makewindows -g Dcitri_chroms_length_tabfile.txt \
-w 10000 -i srcwinnum > Dcitri_chroms_windows.bed

samtools index female_sorted.bam
samtools index male_sorted.bam

for file in $(ls *sorted.bam)
do 
	base=$(basename ${file} "sorted.bam")
	mosdepth -t 4 --no-per-base -b Dcitri_chroms_windows.bed \
	${base} ${base}sorted.bam 
done

echo "moving outputs"
mv *regions.bed.gz /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"