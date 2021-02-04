# Gene Ontology enrichment for differentially expressed and methylated genes

Using the GO annotation file made in the genome paper, it appears there are only 5861 annotated GO terms.

This means only 338/12420 genes used in the differential expression analysis have GO terms and 10/1250 differentially expressed genes have GO terms. This isn't good enough.

I need to re-annotate the genes/transcripts with GO terms. CrowdGO seems the best option for this but there is an issue with dependencies on the qm. I have therefore tried to use a different HPC (University of Leicester, ALICE) and snakemake is pointing to the wrong version of python! I need to speak with the administrators to fix this. I will come back to the GO term analysis later on. Worst case scenario if I contact Marteen (the maker of CrowdGO) maybe he can help me out.

