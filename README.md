# Asian Psyllid Sex-Specific DNA Methylation

This project is a collaboration with Dr. Xiudao Yu and Dr. Zhanjun Lu from Gannan Normal Universtiy, Jiangxi, China. We plan to identifiy genome-wide differences in DNA methylation between male and female *Diaphorina citri*, including the relationship of DNA methylation with gene expression and the DNA methylation differences in the sex chromosomes. This [paper](https://doi.org/10.1111/mec.15216) in aphids shows roughly what we want to achieve.

There is a chromosome level genome with annotation and GO terms from this [paper](https://www.biorxiv.org/content/10.1101/869685v1).

---

## Task 1: Determine the sex chromosome system
Summary: there is an X0 system with chromosome 08 being the X chromosome. See [here](./Identification_Sex_Chromosomes/identification_sex_chromosomes.md) for more details.

---

## Task 2: annotate TEs for later methylation analyses
Summary: TE's successfully annotated using the EDTA pipeline, see [here](./TE_annotation/TE_annotation.md) for more details.

---

## Task 3: differential gene expression between sexes

### To-do list:
- GO enrichment **Need a new GO term annotation**
- Question: do we want to look at selection on these genes? Andrew?

**Current progress:**

For up to date progress see [here](./Differential_expression/Differential_expression.md). 

---

## Task 4: differential DNA methylation between sexes

### To-do list: 
- Genome-wide DNA methylation differences between males and females
    - Overall levels including non-CpG (boxplots)
    - Chromosome level differences (including sex chromosomes) (I like Mather's et al. density plot better than our current one)
    - Mean feature levels (e.g. promotors/exons etc.)
    - Counts in high/medium/low/no bins per feature
    - Exon/intron average levels per gene, **NEED TO DO THIS BEFORE CAN MOVE ON SO CAN FIGURE OUT IF NEED TO SPLIT EXONS INTO BINS, E.G. EXONS 1-3 ... this is proving tricky ...)**
- Differential DNA methylation on the CpG and feature level
    - Location of differentially methylated CpGs and count (could also include a manhatten plot)
    - Relative position in gene (if enriched there)
    - MA / scatter plot 
    - Gene ontology (GO) enrichment of differentially methylated genes 

**Current progress:**

For up to date progress see [here](./Differential_methylation/Differential_methylation.md).

---

## Future tasks:

- Genome-wide relationship between DNA methylation and expression
    - Gene level scatter plots
    - Mean FPKM per binned methylation (line graph)
    - FPKM per highly binned methylation (violin plot)
- Relationship of differentially methylated genes and differentially expressed
    - General scatter of overexpressed genes and hypermethtylated genes etc.
    - Male/female biased methylated/expressed enriched on the X?
    - Overall methylation levels of differentially expressed
    - Overall expression levels of differentially methylated