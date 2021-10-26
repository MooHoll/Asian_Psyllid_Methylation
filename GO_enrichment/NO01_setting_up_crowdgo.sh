# ON ALICE in a single directory: (e.g. bin)
git clone https://gitlab.com/mreijnders/CrowdGO.git

module load python/gcc/35
virtualenv --no-site-packages crowdgo_env
source /scratch/monoallelic/hm257/daphnia/bin/crowdgo_env/bin/activate
python3 -m pip install scipy
python3 -m pip install pandas
python3 -m pip install xgboost
python3 -m pip install imbalanced-learn
python3 -m pip install scikit-learn

module load interproscan/5.46.81
module load snakemake/4.3.1
module load hmmer/3.1b2

wget http://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

gunzip uniprot_sprot.fasta.gz
gunzip uniprot_trembl.fasta.gz
cat uniprot_sprot.fasta >> uniprot_trembl.fasta
mv uniprot_trembl.fasta knowledgebase.fasta
/scratch/monoallelic/hm257/daphnia/bin/diamond makedb -d knowledgebase --in knowledgebase.fasta

# I copied the Pfam-A.hmm file from Edinburgh's qm (Andrés already had it)
hmmpress Pfam-A.hmm

cd tools/cath-tools-genomescan/data
wget http://download.cathdb.info/cath/releases/all-releases/v4_2_0/sequence-data/funfam-hmm3-v4_2_0.lib.gz
gunzip funfam-hmm3-v4_2_0.lib.gz

# Copied goaUniprot_IEA.tab.gz file from Maarten's google drive to my computer then to ALICE
gunzip goaUniprot_IEA.tab.gz

# run this to make sure deepgoplus and the diamond version are compatible (it works!)
cd CrowdGO/tools/deepgoplus/data
/scratch/monoallelic/hm257/daphnia/bin/diamond/diamond makedb --in training_data.fasta --db train_data.dmnd

# TO DO: test run the snake make file with the test dataset

# Edit config.yaml in the CrowdGO directory to the following (all paths relative to the snakemake directory,
# this: /scratch/monoallelic/hm257/daphnia/bin/CrowdGO)

fastaFile: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/example_input/testdata_snakemake.fasta' ## Path to protein fasta se$
ModelFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/models/CrowdGOFull' ## Path to the model used for predicting fast$
outputFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/example_output' ## Path to the Output folder where all predictio$
tmpFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/temp' ## Path to a folder where temporary files can be written
knowledgebaseDBPath: '/scratch/monoallelic/hm257/daphnia/bin/knowledgebase.dmnd' ## Path to the UniProt diam$
PfamDBPath: '/scratch/monoallelic/hm257/daphnia/bin/Pfam-A.hmm' ## Path to the Pfam-A.hmm hmmscan database

mode: 'Full' ## Full, Full-SwissProt, Light, Light-SwissProt

CrowdGOFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO' ## Path where CrowdGO is installed
DeepGOPlusFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/tools/deepgoplus/' ## Path where DeepGOPlus is installed
FunFamsFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/tools/cath-tools-genomescan/' ## Path where FunFams is installe$
InterProScanFolder: 'interproscan.sh' ## Path where InterProScan is installed
Wei2GOFolder: '/scratch/monoallelic/hm257/daphnia/bin/CrowdGO/tools/wei2go/' ## Path where Wei2GO is installed

DiamondBinaryPath: '/scratch/monoallelic/hm257/daphnia/bin/diamond/diamond' ## If the DIAMOND binary is not transfered to a pl$
