# Differential expression between male and female Asian Citrus Psyllid

NOTE: I downloaded the 'clean' data from Novagene rather than the 'raw' data.

## Quality checking of data
**Methods:** Data were quality checked using fastqc v0.11.8 ([Andrews 2010](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). Reads were quality trimmed using CutAdapt v1.18 ([Martin 2011](10.14806/ej.17.1.200)).

**Result:** Need to trim first 10bp for all samples, no big deal. I'm a little concerned that the GC content per read as it shows clear peaks which are consistent across all samples. I've asked the lab to see what they think and they think that it could be contamination, in which case it shouldn't map to the genome and we should be fine. This will only be an issue if the mapping rate is particuarly low, we will see. All fastqc reports of raw data can be found [here](https://www.dropbox.com/sh/ugiz2ipnwt3cm3z/AAB10t28_lZmHDsPYLRDkj-Da?dl=0). All fastqc resports for post-trimmed data can be found [here](https://www.dropbox.com/sh/jdivops93qmp6x5/AABT0c5VJH1N0dX1CzCFK2Yba?dl=0).


---

## Alignment to the reference genome
**Methods:** Samples were aligned to the reference genome Diaci v1.1 ([Hosmani et al. 2020](https://www.biorxiv.org/content/10.1101/869685v1)) using STAR v2.7.3a ([Dobin and Gingeras 2016](10.1002/0471250953.bi1114s51.Mapping)).

**Result:**

---