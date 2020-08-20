#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/mealybugs/analyses/hollie/logs/downloading_data_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`
#---------------------------------------------
# Possibly after 24 hours you have to get a new link as it suddenly goes to "forbidden"
# Also, the samples download with full link name so need renaming but otherwise it's fine

# Sequencing info
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/Report-X101SC20073488-Z02-J001-B1-21.zip?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893593&Signature=hSVbG%2FrOWZdLw0tua6IXR4xVnJs%3D"

# Raw data
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/F_FDSW202227352-1r/F_FDSW202227352-1r_1.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597977528&Signature=ZHLEY%2Bznl08osaN6BpA1pT59L6I%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/F_FDSW202227352-1r/F_FDSW202227352-1r_2.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893637&Signature=etOgWq27wnEXsHcP7S7VRu6vPhI%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/F_FDSW202227352-1r/MD5_F_FDSW202227352-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893637&Signature=L15UuYy4dIhpJAQZibR5oqH1Bl4%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/M_FDSW202227353-1r/MD5_M_FDSW202227353-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893637&Signature=zyeGVb%2BQfI%2BOLQAgfTNGUjMHKGE%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/M_FDSW202227353-1r/M_FDSW202227353-1r_1.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893637&Signature=jHUMHkIzVRd%2F073CjB%2FawAqPl0g%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/1.rawdata/M_FDSW202227353-1r/M_FDSW202227353-1r_2.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893637&Signature=9vboqDaoDpszqzoJw89BJVt4Jp0%3D"

# Clean data
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/F_FDSW202227352-1r/F_FDSW202227352-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=Ou%2Fu%2FV9iZYdsbUMPJoJndN74rWk%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/F_FDSW202227352-1r/F_FDSW202227352-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=hgivrhU0szrXoZ%2FBeND1WK8%2F290%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/F_FDSW202227352-1r/MD5_F_FDSW202227352-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=JzFRFrsy3rgsYCl1PixEy6htnM0%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/M_FDSW202227353-1r/MD5_M_FDSW202227353-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=4KKnZtFdyg3khxA73z3nEc6NKSc%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/M_FDSW202227353-1r/M_FDSW202227353-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=oBlMXyfY55kLbDhEGbNaayaCm5Q%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2019011200887/H101SC20073488/KY_nuohe_JK/X101SC20073488-Z02/X101SC20073488-Z02-J001/2.cleandata/M_FDSW202227353-1r/M_FDSW202227353-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1597893714&Signature=t0gNo8DLkR4wJdsD0MoFjlAurLY%3D"

echo "moving outputs"
mv ./* /data/ross/sequencing/raw/asian_psyllid

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"