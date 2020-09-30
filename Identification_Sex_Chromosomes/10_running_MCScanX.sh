# Making the correct GFF files for chromosome synteny analysis with MCScanX

# The following: chr.list and changing_scaffolds_to_chr.pl are the property of :
# https://github.com/lyy005/Psyllid_chromosome_assembly
#---------------------------------------------------------------
# Copy the below contents to a file called chr.list
#---------------------------------------------------------------
pv1 ScZCZ4B_818 chr1 52474386
pv2 ScZCZ4B_39887 chr2 44032699
pv3 ScZCZ4B_39812 chr3 42045613
pv4 ScZCZ4B_141 chr4 40426472
pvX ScZCZ4B_2870 chr5 38917430
pv5 ScZCZ4B_739 chr6 36706825
pv6 ScZCZ4B_39767 chr7 34277301
pv7 ScZCZ4B_468 chr8 33723484
pv8 ScZCZ4B_39793 chr9 30319696
pv9 ScZCZ4B_16 chr10 22185946
pv10 ScZCZ4B_61 chr11 21967884
pv11 ScZCZ4B_504 chr12 18660665
pv12 ScZCZ4B_39851 chr13 2654893
ap1 Scaffold_20849 Chr1 0 170740645 0,115,194
ap2 Scaffold_21967 Chr2 0 119541763 134,134,134
ap3 Scaffold_21646 Chr3 0 42333646 205,83,76
apX Scaffold_21773 ChrX 0 132544852 239,192,0
rm1 NC_040877.1
rm2 NC_040878.1
rm3 NC_040880.1
rmX NC_040879.1
#chr - Scaffold_20960 Sf1 0 1832872 black
#chr - Scaffold_21109 Sf2 0 1247914 black

#---------------------------------------------------------------
# Copy below contents to script called changing_scaffolds_to_chr.pl
#---------------------------------------------------------------
#!/usr/bin/perl -w
use strict;

die "usage: perl $0 [ old gff file ] [ chr list ] [ new gff file ]\n" unless (@ARGV == 3);

open GFF,$ARGV[0] or die "$!\n";
open LST,$ARGV[1] or die "$!\n";
open OUT,">$ARGV[2]" or die "$!\n";

my %hash;
print "#Scaffold_ID\tchr_ID\n";
while (my $line=<LST>) {
        next if ($line =~ /^\#/);
        chomp $line;
        my @line = split(/\s+/,$line);
        $hash{$line[1]} = $line[0];
print "$line[1]\t$line[0]\n";
}
close LST;
#---------------------------------------------------------------

perl changing_scaffolds_to_chr.pl Pachypsylla_venusta.gff chr.list Pachypsylla_venusta_new.gff3

#---------------------------------------------------------------

grep "gene" Pachypsylla_venusta_new.gff > look1
cut -f1,4,5,9 look1 > look2
sed -i 's/ID=//g' look2 
sed -i 's/;//g' look2 
awk '{ print $1 " " $4 " " $2 " " $3 }' look2 > Pachypsylla_venusta_edited.gff

grep "gene" Dcitr_OGSv3.0_beta.gff3 > look1
cut -f1,4,5,9 look1 > look2
sed -i 's/3.0sc//g' look2
sed -i 's/ID=//g' look2 
sed -i 's/;.*//g' look2
awk '{ print $1 " " $4 " " $2 " " $3 }' look2 > Dcitri_edited.gff

cat Pachypsylla_venusta_edited.gff Dcitri_edited.gff > Dcitri_to_Pvenusta.gff

#cat Dcitri_to_Pvenusta_test.gff | sort | uniq > Dcitri_to_Pvenusta_singlecopy.gff

#---------------------------------------------------------------



#### TIME TO STOP:
# Problem:
# The gene name apparently does not correspond in both the .gff and the blast result
# this is why I am getting this: 0 homologous pairs imported (71778 discarded)
# See https://www.programmersought.com/article/1848716053/
# Ahh ok I've found it, in the GFF the second column should be the gene name not the word gene duh!