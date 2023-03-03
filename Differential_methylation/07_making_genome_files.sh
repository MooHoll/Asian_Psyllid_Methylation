############################################################################
### Making genome information files
############################################################################

# Use one genome-wide cytosine report from methylation_extraction from bismark
# to create a file which consists of 
# just the scaffold name and the CpG position in a text file
# Only take + stranf coorsinate otherwise get each CpG counted twice
grep "+" F1-1_BDLM200001629-1A_1.clean_bismark_bt2.deduplicated.CpG_report.txt > plus_only.txt
cut -f1,2 plus_only.txt > total_cpgs_in_genome.txt

#---------------------------------------------------

# Make gene name and start and end positions
grep "gene" Dcitr_OGSv3.0_beta.gff3  > genes.txt
cut -f1,4,5,7,9 genes.txt > genes_next.txt
sed 's/ID=//g' genes_next.txt > new.txt
sed 's/;Name=.*$//g' new.txt > new1.txt
echo -e "chr\tstart\tend\tstrand\tgene_id" | cat - new1.txt > genes_with_start_and_end.txt









