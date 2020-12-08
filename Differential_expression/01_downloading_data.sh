#$ -V
#$ -cwd
#$ -j y
#$ -o /data/ross/misc/analyses/asian_psyllid/logs/downloading_data_$JOB_ID.o

set -e

SCRATCH=/scratch/$USER
mkdir -p $SCRATCH
cd $SCRATCH

start=`date +%s`
#---------------------------------------------

#RNA-Seq female
#wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-1_FRAS202200995-1r/F2-1_FRAS202200995-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503367&Signature=B2ttIV4RGoqrCuYxOjbOmgn9%2BUo%3D"
#wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-1_FRAS202200995-1r/F2-1_FRAS202200995-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503553&Signature=m8ZKAWl4ckeEETUJURqeAT8ibfY%3D"
#wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-2_FRAS202200996-1r/F2-2_FRAS202200996-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503795&Signature=%2FpvIi0npp3uiqOn%2B%2BxYLbfGB664%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-2_FRAS202200996-1r/F2-2_FRAS202200996-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503842&Signature=PYP%2BRG%2Be8EhiC2GcroLmLe%2Fd%2BzY%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-3_FRAS202200997-1r/F2-3_FRAS202200997-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508520&Signature=tKYm%2Fn000Q0zdV72sNe55eP6Zew%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-3_FRAS202200997-1r/F2-3_FRAS202200997-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508549&Signature=8W5BE2uKYAFgWSh8Nu%2FDNcOQwsY%3D"

#md5s for RNA-seq female
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-1_FRAS202200995-1r/MD5_F2-1_FRAS202200995-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503581&Signature=YtR8LIABWkQhCsugSPGTxwWp3OI%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-2_FRAS202200996-1r/MD5_F2-2_FRAS202200996-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607503862&Signature=BJcn4YYrLSDzGQB2k4EC8jwYwqg%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/F2-3_FRAS202200997-1r/MD5_F2-3_FRAS202200997-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508606&Signature=xBEU2nsJPNKfe2zaP3hYKhk1srA%3D"

#RNA-Seq male
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-1_FRAS202200998-1r/M2-1_FRAS202200998-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508667&Signature=AAoOC91R83EOkKzmPt0Pc4ZlDiA%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-1_FRAS202200998-1r/M2-1_FRAS202200998-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508675&Signature=VrDHkL68%2FRI50z%2BkUDZ7hUshF9o%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-2_FRAS202200999-1r/M2-2_FRAS202200999-1r_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508725&Signature=2QzpsAa%2FE6kwe1EYqY47PT%2FtlgM%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-2_FRAS202200999-1r/M2-2_FRAS202200999-1r_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508732&Signature=yq%2BTwWurqGazyWmsnA9bGk3SUHQ%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-3_FRAS202201000-1b/M2-3_FRAS202201000-1b_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508826&Signature=mXShHKtHx4sku2n59LxfjLehegA%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-3_FRAS202201000-1b/M2-3_FRAS202201000-1b_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508835&Signature=7gtPEfWTCbqj4Dtvo5gsQT1FhWg%3D"

#md5s for RNA-seq male
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-1_FRAS202200998-1r/MD5_M2-1_FRAS202200998-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508685&Signature=BD6vBrk%2FmGYrdlL1U6cCO7k6Oqk%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-2_FRAS202200999-1r/MD5_M2-2_FRAS202200999-1r.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508746&Signature=iVe1vvxDmrHxHDkIyb%2Bcux5I0ks%3D"
wget "http://novo-intelligent-tj.oss-cn-beijing.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS3300/X101SC20101161-Z01/X101SC20101161-Z01-J002/2.cleandata/M2-3_FRAS202201000-1b/MD5_M2-3_FRAS202201000-1b.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508844&Signature=0iqYfDITZNQZIK5Fx7IQgeQfeqQ%3D"

#WGBS
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-1_BDLM200001629-1A/F1-1_BDLM200001629-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508954&Signature=xLToILswsNntClbI7E%2B8ckCsRFE%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-1_BDLM200001629-1A/F1-1_BDLM200001629-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607508981&Signature=XS%2FPvmHqTf062tRDXrMne44yM7g%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-2_BDLM200001630-1A/F1-2_BDLM200001630-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509144&Signature=jJ8eDEj4%2FjzFGZZ9gWLw7X8KVlY%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-2_BDLM200001630-1A/F1-2_BDLM200001630-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509158&Signature=C5hYdq%2BX6OBNgO1CPL%2Fv7mRTVRQ%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-3_BDLM200001631-1A/F1-3_BDLM200001631-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509247&Signature=7cGFQkmUmwEiRVyD8A%2BfNXse9b0%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-3_BDLM200001631-1A/F1-3_BDLM200001631-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509258&Signature=KMCLpONpoVJafi8mxZZUY674skk%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-1_BDLM200001632-1A/M1-1_BDLM200001632-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509363&Signature=%2FXfXXytwsBjFkgWwhYL1xmKS%2FP8%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-1_BDLM200001632-1A/M1-1_BDLM200001632-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509378&Signature=R435ThynhK9fr2kr7i6TZXv6m68%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-2_BDLM200001633-1A/M1-2_BDLM200001633-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509403&Signature=Bs8PveuHjlwHx420Wv7N3XhTX1Y%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-2_BDLM200001633-1A/M1-2_BDLM200001633-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509412&Signature=GyF1VrLfSxRPZtocZuVsz940xHM%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-3_BDLM200001634-1A/M1-3_BDLM200001634-1A_1.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509456&Signature=I3KbaN2FD5uSAb2GE%2BSrlgOKllY%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-3_BDLM200001634-1A/M1-3_BDLM200001634-1A_2.clean.fq.gz?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509466&Signature=aLu%2Btd6T6NVIDArphRtk3CS%2FXBY%3D"

#MD5s
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-1_BDLM200001629-1A/MD5_F1-1_BDLM200001629-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509005&Signature=9HzzOV4%2Bld5OdzZVsDpdkCUR8Xg%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-2_BDLM200001630-1A/MD5_F1-2_BDLM200001630-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509201&Signature=LMsoo%2B212EtyYyPPQAoxthJVaUI%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/F1-3_BDLM200001631-1A/MD5_F1-3_BDLM200001631-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509284&Signature=r5VK%2FqeloaUpgeVKJJNEfwXXMxY%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-1_BDLM200001632-1A/MD5_M1-1_BDLM200001632-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509387&Signature=sil7I5p3UbqniAqYKs0mPsJgCnk%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-2_BDLM200001633-1A/MD5_M1-2_BDLM200001633-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509422&Signature=2qXbMeUo3icKJZqb4K50tcns%2B6s%3D"
wget "http://novo-intelligent-nj.oss-cn-hangzhou.aliyuncs.com/css_home/CP2020082000115/H101SC20101161/RSCS5800/X101SC20101161-Z02/X101SC20101161-Z02-J001/2.cleandata/M1-3_BDLM200001634-1A/MD5_M1-3_BDLM200001634-1A.txt?OSSAccessKeyId=LTAIBg5ZREoSx2z9&Expires=1607509476&Signature=VRT65AjkPaV9vhaLlcKlhvE1Jpw%3D"

#---------------------------------------------

echo "moving outputs"
mv ./* /data/ross/sequencing/raw/asian_psyllid/raw_data

echo "a clean directory is a happy directory"
rm -r $SCRATCH/*

#---------------------------------------------

end=`date +%s`
runtime=$((((end-start)/60)/60))
echo "runtime:"$runtime"hrs"