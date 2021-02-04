# Making a bed-type file of the annotation for manipulation

# keep only the longest transcript
conda activate hmarshall
conda install -c bioconda agat
agat_sp_keep_longest_isoform.pl -gff Dcitr_OGSv3.0_beta.gff3 -o Dcitr_OGSv3.0_beta_longestIsoformOnly.gff3
agat_sp_add_introns.pl -gff Dcitr_OGSv3.0_beta_longestIsoformOnly.gff3 -o Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons.gff3

# make file to get intergenic regions
grep -v "#" Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons.gff3 > annotation1
sed 's/ID=.*Parent=//g' annotation1 > annotation2
sed 's/ID=//g' annotation2 > annotation3
sed 's/;Name=.*$//g' annotation3 > annotation4
sed 's/\.1.*$/\.1/1' annotation4 > Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons.txt

# Fix the TE annotation file as well
sed 's/;Sequence_ontology.*$//g' Diaci_v3.0.ref.fa.mod.EDTA.intact.gff3 > te_annot1
sed 's/ID=.*Classification=//g' te_annot1 > Diaci_v3.0.ref.fa.mod.EDTA.intact.txt


