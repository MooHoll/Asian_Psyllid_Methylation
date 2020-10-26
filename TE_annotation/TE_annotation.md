# Annotate TEs in the Asian Citrus Psyllid

We want to annotate TEs in *D. citri* so we can check the DNA methylation levels.

There are two ways I can do this: 

1. follow the RepeatMaskers:RepeatModeler pipeline [here](https://github.com/RossLab/Sex-Specific_Methylation_P.citri/tree/master/Transposable_elements)

2. Try a new pipeline which integrates this pipeline with others for a comprehensive identification of TEs, it also categorises the TEs which looks cool, paper [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y) and the software [here](https://github.com/oushujun/EDTA).

---

## Running EDTA for TE annotation

### Create EDTA environment

Install as shown in the softwar github [here](https://github.com/oushujun/EDTA).

    conda create -n EDTA_env
    conda activate EDTA_env
    conda config --env --add channels anaconda --add channels conda-forge --add channels bioconda
    conda install -n EDTA_env -y cd-hit repeatmodeler muscle mdust blast openjdk perl perl-text-soundex multiprocess regex tensorflow=1.14.0 keras=2.2.4 scikit-learn=0.19.0 biopython pandas glob2 python=3.6 tesorter genericrepeatfinder genometools-genometools ltr_retriever ltr_finder numpy=1.16.4
    git clone https://github.com/oushujun/EDTA

### Make a genes bed file for *D.citri*

    grep "gene" Dcitr_OGSv3.0_beta.gff3 > genes
    cut -f 4,5,9 genes > genes_cut
    sed -i 's/ID=//g' genes_cut 
    sed -i 's/;.*//g' genes_cut 
    awk '{ print $3 "\t" $1 "\t" $2}' genes_cut > Dcitri_genes.bed

### Output files

Everything ran well with no errors or warnings. There are a LOT of output files, here are links to the main ones of interest: 

[Diaci_v3.0.ref.fa.mod.EDTA.TElib.fa](./files/Diaci_v3.0.ref.fa.mod.EDTA.TElib.fa) - fasta seq of the different TEs found <br/>
[Diaci_v3.0.ref.fa.mod.EDTA.intact.gff3](./files/Diaci_v3.0.ref.fa.mod.EDTA.intact.gff3) - gff file of the TEs found <br/>
[Diaci_v3.0.ref.fa.mod.EDTA.TEanno.gff3](./files/Diaci_v3.0.ref.fa.mod.EDTA.TEanno.gff3) - gff file of the TEs found including fragmented TEs