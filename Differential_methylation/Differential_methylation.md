# Differential DNA methylation between male and female Asian Citrus Psyllid


## Quality checking of data
**Methods:** Data were quality checked using fastqc v0.11.8 ([Andrews 2010](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).

**Results:**
The data look beatutiful, no need for any trimming (already done by novagene I believe as reads are 140bp). Fastqc reports can be found [here](https://www.dropbox.com/sh/x77ogcoe0um250i/AAB-PO2nCbTiKzgdi9fPB8aEa?dl=0).

---

## Alignment to the reference genome
**Methods:** Samples were aligned to the reference genome Diaci v3.0 ([Hosmani et al. 2020](https://www.biorxiv.org/content/10.1101/869685v1)) and also the lambda genome to determine the bisulfite conversion rate using Bismark v0.20.0 ([Krueger and Andrews 2011](https://pubmed.ncbi.nlm.nih.gov/21493656/)).

**Results:**
0% reads have aligned to the reference genome or the lambda genome in paired-end mode. Xiudao sent the methods used by Novogene and it all looked standard. Novogene suggested aligning using the --dovetail argument and this has done nothing. I also tried a different lambda genome and still nada. Following [this](https://github.com/FelixKrueger/Bismark/blob/master/Docs/FAQ.md#mapping-strategies-for-paired-end-data) I have now successfully aligned the reads WOO! I have used single-end mode and we have about 1% of reads aligning to the lambda genome and the conversion looks good. Using just single-end mode (and --pbat for read 2) we get about 30% alignment rate to Diaci v3.0. By decreasing the alignment stringency slightly the alignment rate has increased to about 50% which is much better.

Alignment reports can be found [here](https://www.dropbox.com/sh/nnbg3ff8mqrde29/AAAzwR98UImg7NNKR_Qz9wgOa?dl=0). 

At a first look it appears *D. citri* has about 1% CpG methylation and 0% in a non-CpG context. This is pretty standard for an insect and a useable amount. Exciting!

---

**Currently running deduplication.**