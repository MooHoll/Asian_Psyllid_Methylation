# ON ALICE in a single directory: (e.g. bin)
git clone https://gitlab.com/mreijnders/CrowdGO.git

module load python/gcc/3.6.4
virtualenv --no-site-packages crowdgo_env
source /scratch/monoallelic/hm257/daphnia/bin/crowdgo_env/bin/activate
pip3 install scipy
pip3 install pandas
pip3 install xgboost
pip3 install imbalanced-learn
pip3 install scikit-learn

module load interproscan/5.46.81
module load snakemake/4.3.1
module load hmmer/3.1b2

wget http://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

gunzip uniprot_sprot.fasta.gz
gunzip uniprot_trembl.fasta.gz # RUN FROM HERE
cat uniprot_sprot.fasta >> uniprot_trembl.fasta
mv uniprot_trembl.fasta knowledgebase.fasta
diamond makedb -d knowledgebase --in knowledgebase.fasta