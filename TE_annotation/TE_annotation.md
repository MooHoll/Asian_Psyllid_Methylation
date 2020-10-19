# Annotate TEs in the Asian Citrus Psyllid

We want to annotate TEs in *D. citri* so we can check the DNA methylation levels.

There are two ways I can do this: 

1. follow the RepeatMaskers:RepeatModeler pipeline [here](https://github.com/RossLab/Sex-Specific_Methylation_P.citri/tree/master/Transposable_elements)

2. Try a new pipeline which integrates this pipeline with others for a comprehensive identification of TEs, it also categorises the TEs which looks cool, paper [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y) and the software [here](https://github.com/oushujun/EDTA)

---

## Create an environment for the EDTA doftware
    conda create -n EDTA_env
    conda activate EDTA_env
    conda config --env --add channels anaconda --add channels conda-forge --add channels bioconda
    conda install -n EDTA_env -y cd-hit repeatmodeler muscle mdust blast openjdk perl perl-text-soundex multiprocess regex tensorflow=1.14.0 keras=2.2.4 scikit-learn=0.19.0 biopython pandas glob2 python=3.6 tesorter genericrepeatfinder genometools-genometools ltr_retriever ltr_finder numpy=1.16.4
    git clone https://github.com/oushujun/EDTA
    ./EDTA/EDTA.pl