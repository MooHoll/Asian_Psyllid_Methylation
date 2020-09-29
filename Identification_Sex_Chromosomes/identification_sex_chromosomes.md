# Determine the sex chromosome system

There is a chromosome level assembly of *D. citri* but the sex chromosomes have not been identified. Xiudao has therefore carried out whole genome resequencing of male and female *D. citri* so we can attempt to determine whether the system is XO or XY.

Genome and annotation: the chromosome level assembly paper is avaliable [here](https://www.biorxiv.org/content/10.1101/869685v1). As it is still a pre-print and not published ther authors have not yet released the full version. I have emailed the authors and obtained version 3 plus a current annotation file.

I will follow Kamil's coverage pipeline to see if we have clear coverage differences in one particular chromosome set in males vs females: [https://github.com/RossLab/PGE/blob/master/sex_chromosome_p_citry.md.](https://github.com/RossLab/PGE/blob/master/sex_chromosome_p_citry.md).

---

Following Kamil's pipeline above I have successfully identified that the X chromosome is chromosome 8. I have now also ran a windowed coverage approach to look for the Y chromosome. I have found no clear distributions >=1 male2femaleLogRatio meaning there is no divereged Y present. The system is most likely X0 however, I cannot rule out a newly evolved Y chromosome using this method. In addition to the figures generated below I also took all windows throughout the genome which had 0 / very low coverage in females and plotted the coverage for males for these windows to manually check if there are any small regions not visable on the graphs and there are none. This means there is truely no divergant Y chromosome.

<img src="./images/coverage_plot.jpeg" height="240">
Fig.1: Whole chromosome coverage analysis shows chromosome 8 is the X chromosome.<br/>

<img src="./images/all_chr_distribution.jpeg" height="240"> <img src="./images/per_chr_distribution.jpeg" height="240"> <br/>
Fig.2: coverage analysis of each chromosome split into 10,000bp windows shows no distributions are present above the male to female logRatio 1.0 meaning there are no regions which are present in the male only, suggesting there is no Y chromosome present.<br/>

<img src="./images/08.png" height="240">
Fig.3: A closer look at chromosome 8 shows a distribution again with no peaks >1.0 meaning there are no male only windows. <br/>

**Question: is this analysis good enough to confirm X0 or is there more we can do, with the RNA-Seq for example?**

**Answer: We think so, the original genome was based mostly on male data and we should therefore see a clear Y if there was one.**

---

## Checking synteny of the X chromosome of D. citri with another psyllid

The psyllid *Pachypsylla venusta* has recently been sequneced to chromosome level and the X chromosome identified, see the publication [here](https://academic.oup.com/mbe/article/37/8/2357/5820017). I will check for chromosome synteny between *D. citri* and *P. venusta* to check how similar the X chromosomes are.

The genome can be found [here](https://www.ncbi.nlm.nih.gov/genome/?term=txid38123[orgn]), version: Pven_dovetail. The annotation files and protein fasta files can be found [here](https://github.com/lyy005/Psyllid_chromosome_assembly/tree/master/step0_genome_annotation_files). 


**NOTES:** Currently running the blast, next need to make the .gff file into the correct format for MCScanX, see [here](http://chibba.pgml.uga.edu/mcscan2/documentation/manual.pdf).
