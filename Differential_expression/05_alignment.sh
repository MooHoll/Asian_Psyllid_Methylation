#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/alignment_RNA_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`

#---------------------------------------------

echo "copying genome data in"
rsync /data/ross/misc/analyses/asian_psyllid/genome/Diaci_v3.0.ref.fa ./
rsync /data/ross/misc/analyses/asian_psyllid/genome/Dcitr_OGSv3.0_beta.gff3 ./

mkdir genome_dir

echo "prepare the genome"
    STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles Diaci_v3.0.ref.fa \
    --sjdbGTFfile Dcitr_OGSv3.0_beta.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --runThreadN 8 \
    --genomeDir ./genome_dir 

echo "copying RNA-seq in"
rsync /data/ross/misc/analyses/asian_psyllid/raw_data/trimmed_RNA-Seq/*fq.gz ./

echo "align to genome"
for file in $(ls *_1.fq.gz)
do
    base=$(basename $file "_1.fq.gz")
    STAR \
    --runThreadN 32 \
    --genomeDir ./genome_dir \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}_1.fq.gz ${base}_2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}
done

echo "moving outputs"
mv ./* /data/ross/misc/analyses/asian_psyllid/3_expression

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$(((end-start)/60))
echo "runtime:"$runtime"mins"
