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
-w 10000 -i srcwinnum > Dcitri_chroms_windows.txt

#awk '{ print $4 " " $2 " " $3}' Dcitri_chroms_windows.txt > Dcitri_chroms_windows.bed

for file in $(ls *sorted.bam)
do 
	base=$(basename ${file} "sorted.bam")
	samtools depth -b Dcitri_chroms_windows.txt ${file} | \
    awk '/BEGIN/{scf='DC3.0sc00'; coverage_sum = 0; }{ if( scf != $1 ){ print scf "\t" coverage_sum; scf = $1; coverage_sum = $3 } else { scf = $1; coverage_sum += $3} }' > ${base}window_depth.txt
done

echo "moving outputs"
mv *depth.txt /data/ross/misc/analyses/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
