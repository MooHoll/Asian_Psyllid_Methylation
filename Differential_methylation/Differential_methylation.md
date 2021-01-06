# Differential DNA methylation between male and female Asian Citrus Psyllid


## Quality checking of data
**Methods:** Data were quality checked using fastqc v0.11.8 ([Andrews 2010](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).

**Results:**
The data look beatutiful, no need for any trimming (already done by novagene I believe as reads are 140bp). Fastqc reports can be found [here](https://www.dropbox.com/sh/x77ogcoe0um250i/AAB-PO2nCbTiKzgdi9fPB8aEa?dl=0).

---

## Alignment to the reference genome
**Methods:** Samples were aligned to the reference genome Diaci v3.0 ([Hosmani et al. 2020](https://www.biorxiv.org/content/10.1101/869685v1)) and also the lambda genome to determine the bisulfite conversion rate using Bismark v0.20.0 ([Krueger and Andrews 2011](https://pubmed.ncbi.nlm.nih.gov/21493656/)).

**Results:**
We have a problem, 0% reads have aligned to the reference genome or the lambda genome, see [here](https://www.dropbox.com/sh/nnbg3ff8mqrde29/AAAzwR98UImg7NNKR_Qz9wgOa?dl=0). 

Xiudao has sent the methods used by Novogene and it all looks standard. Novogene suggested aligning using the --dovetail argument and this has done nothing, I have also tried a different lambda genome, also still nothing.

*I am now running the alignment in single ended mode as there may be a problem with the paired-end nature of the data, see [here](https://github.com/FelixKrueger/Bismark/blob/master/Docs/FAQ.md#mapping-strategies-for-paired-end-data)*