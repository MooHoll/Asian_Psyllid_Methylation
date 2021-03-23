# Identification of DNMT genes to check for differences in expression levels

According to the supplementary of [Bewick et al. (2017)](https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msw264) *D. citri* has two possible DNMT1 genes (XP_008474257.1; XP_008481189.1) and one possible DNMT2 gene (XP_008470452.1) but no DNMT3 gene. These are old annotations from the previous genome so I need to see what they match in the new genome. 

Download protein fasta sequences of the above from NCBI.

Make a Blast database from the protein sequences of the recent genome annotation:

`/ceph/software/blast/ncbi-blast-2.8.0+/bin/makeblastdb 
-in Dcitr_OGSv3.0_beta_pep.fa 
-dbtype prot 
-parse_seqids 
-out D_citri_protein_db`

Blast each of the three proteins

`/ceph/software/blast/ncbi-blast-2.8.0+/bin/blastp -query DNMT2.fa 
-db D_citri_protein_db 
-num_threads 10 
-max_target_seqs 1 
-outfmt 6 
-evalue 1e-3 > DNMT2.txt`

DNMT1-1 : XP_008474257.1 : Dcitr08g10610.1.2
DNMT1-2 : XP_008481189.1 : Dcitr08g05090.1.1
DMMT2 : XP_008470452.1 : Dcitr07g04270.1.1