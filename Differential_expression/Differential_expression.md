# Differential expression and alternative splicing between male and female Asian Citrus Psyllid

NOTE: I downloaded the 'clean' data from Novagene rather than the 'raw' data.

## Quality checking of data
**Methods:** Data were quality checked using fastqc v0.11.8 ([Andrews 2010](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). Reads were quality trimmed using CutAdapt v1.18 ([Martin 2011](10.14806/ej.17.1.200)).

**Result:** Need to trim first 10bp for all samples, no big deal. I'm a little concerned that the GC content per read as it shows clear peaks which are consistent across all samples. I've asked the lab to see what they think and they think that it could be contamination, in which case it shouldn't map to the genome and we should be fine. This will only be an issue if the mapping rate is particuarly low, we will see. All fastqc reports of raw data can be found [here](https://www.dropbox.com/sh/ugiz2ipnwt3cm3z/AAB10t28_lZmHDsPYLRDkj-Da?dl=0). All fastqc resports for post-trimmed data can be found [here](https://www.dropbox.com/sh/jdivops93qmp6x5/AABT0c5VJH1N0dX1CzCFK2Yba?dl=0).

---

## Alignment to the reference genome
**Methods:** Samples were aligned to the reference genome Diaci v3.0 ([Hosmani et al. 2020](https://www.biorxiv.org/content/10.1101/869685v1)) using STAR v2.7.3a ([Dobin and Gingeras 2016](10.1002/0471250953.bi1114s51.Mapping)).

**Results:**
Alignment rate is 50-60% see [here](https://www.dropbox.com/sh/7uhw69wcytcw11w/AACtj_yeM2dUChEr7yHeNZtXa?dl=0). Sequencing company was happy with this, after chatting to others this may be because whole bodies were used which includes the gut etc. and so probably a lot of bacteria were also sequenced. 

---

## Differential expression analysis
**Methods:** RSEM v1.3.3 ([Li and Dewey 2011](http://www.biomedcentral.com/1471-2105/12/323)) was used to calculate gene and transcript abundances and DESEQ2 v1.28.1 ([Love et al. 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)) was used to determine differentially expressed genes between males and females. A gene was considered differentially expressed if the corrected p-value was &lt;0.05 (adjused for multiple testing using the Benjamini-Hochberg procudure ([Benjamini and Hochberg 1995](http://www.jstor.org/stable/2346101%5Cnhttp://about.jstor.org/terms)) and the log2 fold-change was &gt;1.5.


**Results:** See [here](https://www.dropbox.com/sh/mlx50j6znpqbu0b/AAAKgn_sLn7P0t-cftRPMUJMa?dl=0) for the curent graphs and output files.


---

## Differential alternative splicing analysis
**Methods:** DEXSeq vX (REF) implemented by IsoformSwitchAnalyzeR vXX (REF) was used to identify differentially alternativly spliced genes between the sexes. A gene was considered differentially alternativly spliced if the corrected p-value was &lt;0.05 (adjused for multiple testing using the Benjamini-Hochberg procudure ([Benjamini and Hochberg 1995](http://www.jstor.org/stable/2346101%5Cnhttp://about.jstor.org/terms)) and the absolute isoform usage difference was &gt;10%.

**Results:** DOING.