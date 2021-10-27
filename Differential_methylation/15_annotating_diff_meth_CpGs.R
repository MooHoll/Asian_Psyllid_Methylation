## -------------------------------------------------------------------------
# Filtering diff meth CpGs from methylkit with weighted meth of features
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation")
library(sqldf)
library(readr)
library(doBy)
library(ggplot2)
library(dplyr)
library(data.table)

annotation <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files/Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters_plusTEs_plusIntergenic.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

diff_meth_sites <- read_csv("F_vs_M__DMRs_min10percentDiff_qval0.01.csv")
diff_meth_sites <- diff_meth_sites[,c(2,3,8)]
colnames(diff_meth_sites) <- c("chr","cpg_position","methylation_diff")

# Add column for hypermeth sex
diff_meth_sites$hypermethylated <- "female"
diff_meth_sites$hypermethylated[diff_meth_sites$methylation_diff >0] <- "male"

# How many hypermeth in each sex?
nrow(diff_meth_sites[diff_meth_sites$hypermethylated=="female",]) #443
nrow(diff_meth_sites[diff_meth_sites$hypermethylated=="male",]) #320

# Goodness of fit
observed = c(443, 320)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) # X-squared = 19.828, df = 1, p-value = 8.472e-06


## -------------------------------------------------------------------------
# Annotate the differentially methylated CpGs with genomic features

output <- sqldf("SELECT sample.chr,
                    sample.cpg_position,
                    sample.methylation_diff,
                    sample.hypermethylated,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM diff_meth_sites AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg_position >= annot.start AND sample.cpg_position <= annot.end)")

output <- output[,-1]

## -------------------------------------------------------------------------
# Where are these CpGs
## -------------------------------------------------------------------------

not_in_feature <- output[is.na(output$gene_id),] #0 not in feature
output$feature <- as.factor(output$feature)

# Only a couple of splice sites so remove as well
output <- output[!output$feature == "non_canonical_five_prime_splice_site",]
output <- output[!output$feature == "non_canonical_three_prime_splice_site",]

# Remove mRNA and CDS not really informative
output <- output[!output$feature == "CDS",]
output <- output[!output$feature == "mRNA",]

# Remove gene as this is already accounted for by exon/intron etc
output <- output[!output$feature == "gene",]

# Note: 1118 annotations from 763 positions as some positions fall over multiple annotations

ggplot(output, aes(x=feature, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Exon","Intron","TE","Intergenic"),
                   limits =c("intron","exon","intergenic","three_prime_UTR","TE","promoter","five_prime_UTR"))


## -------------------------------------------------------------------------
# Significant CpGs by chromosome
## -------------------------------------------------------------------------

output_labelled <- output
output_labelled$chr[output_labelled$chr=="DC3.0sc08"] <- "X"

ggplot(output_labelled, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))

output_labelled_exons <- output_labelled[output_labelled$feature=="exon",]
output_labelled_introns <- output_labelled[output_labelled$feature=="intron",]

ggplot(output_labelled_exons, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Significant CpGs in Exons")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))

ggplot(output_labelled_introns, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Significant CpGs in Introns")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))

## -------------------------------------------------------------------------
# average number diff CpGs by feature
## -------------------------------------------------------------------------

head(output)
output_count <- summaryBy(cpg_position ~ feature + gene_id , data=output, FUN=length)
output_count <- output_count[!output_count$feature=="intergenic",]

ggplot(output_count, aes(x=feature, y=cpg_position.length))+
  geom_boxplot()+
  xlab("Feature")+
  ylab("Number of Significant CpGs")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 14))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","exon","intron","TE"),
                   labels = c("Promoter","5' UTR","3' UTR","Exon","Intron","TE"),
                   limits =c("intron","exon","three_prime_UTR","TE","promoter","five_prime_UTR"))


# Here we are saying a min of 2 diff meth cpg per feature to count
intron_data <- output[output$feature=="intron",]
number_diff_cpgs_per_intron <- dplyr::count(intron_data, gene_id)
mean(number_diff_cpgs_per_intron$n) #1-7, mean = 1.349, median = 1
hist(number_diff_cpgs_per_intron$n)
nrow(number_diff_cpgs_per_intron[number_diff_cpgs_per_intron$n >=2,]) #55
introns_with_2_cpgs <- subset(number_diff_cpgs_per_intron, n >=2)

exon_data <- output[output$feature=="exon",]
number_diff_cpgs_per_exon <- dplyr::count(exon_data, gene_id)
median(number_diff_cpgs_per_exon$n) #1-7, mean = 1.256, median = 1
hist(number_diff_cpgs_per_exon$n)
nrow(number_diff_cpgs_per_exon[number_diff_cpgs_per_exon$n >=2,]) #40
exons_with_2_cpgs <- subset(number_diff_cpgs_per_exon, n >=2)

UTR3_data <- output[output$feature=="three_prime_UTR",]
number_diff_cpgs_per_UTR3 <- dplyr::count(UTR3_data, gene_id)
range(number_diff_cpgs_per_UTR3$n) #1-6, mean = 1.5, median = 1
hist(number_diff_cpgs_per_UTR3$n)
nrow(number_diff_cpgs_per_UTR3[number_diff_cpgs_per_UTR3$n >=2,]) #7
UTR3_with_2_cpgs <- subset(number_diff_cpgs_per_UTR3, n >=2)

UTR5_data <- output[output$feature=="five_prime_UTR",]
number_diff_cpgs_per_UTR5 <- dplyr::count(UTR5_data, gene_id)
median(number_diff_cpgs_per_UTR5$n) #1-2, mean = 1.13, median = 1
hist(number_diff_cpgs_per_UTR5$n)
nrow(number_diff_cpgs_per_UTR5[number_diff_cpgs_per_UTR5$n >=2,]) #2 
UTR5_with_2_cpgs <- subset(number_diff_cpgs_per_UTR5, n >=2)

promoter_data <- output[output$feature=="promoter",]
number_diff_cpgs_per_prom <- dplyr::count(promoter_data, gene_id)
range(number_diff_cpgs_per_prom$n) #1-2, mean = 1.21, median = 1
hist(number_diff_cpgs_per_prom$n)
nrow(number_diff_cpgs_per_prom[number_diff_cpgs_per_prom$n >=2,]) #5
proms_with_2_cpgs <- subset(number_diff_cpgs_per_prom, n>=2)

# Conclusion focus only on the exon and intron genes

## -------------------------------------------------------------------------
# Take a look at the TEs before we dive into genes (exons/introns)
## -------------------------------------------------------------------------
TEs <- subset(output, feature =="TE")
head(TEs)
TEs$uniquie_identified <- paste0(TEs$start, TEs$gene_id)
length(unique(TEs$uniquie_identified)) #23 TEs

# Filter out all TEs overlapping a significant promotor or exon
# First make a list of CpGs which are diff meth in TEs
cpgs_in_TEs <- unique(TEs$cpg_position)
cpgs_in_introns <- unique(output$cpg_position[output$feature=="exon"])
cpgs_in_exons <- unique(output$cpg_position[output$feature=="intron"])

cpgs_in_both <- cpgs_in_exons[cpgs_in_exons %in% cpgs_in_introns] #28 cpgs in both introns and exons ... urgh

cpgs_unique_to_TEs <- cpgs_in_TEs[!cpgs_in_TEs %in% cpgs_in_introns]
cpgs_unique_to_TEs <- cpgs_unique_to_TEs[!cpgs_in_TEs %in% cpgs_in_exons] # only 14

TEs_subset <- TEs[TEs$cpg_position %in% cpgs_unique_to_TEs,]

# Filter by cpgs per TE
unique_TEs_only <- unique(TEs_subset$uniquie_identified)
number_diff_cpgs_per_TE<- dplyr::count(TEs_subset,uniquie_identified )
TEs_with_2_cpgs <- subset(number_diff_cpgs_per_TE, n >=2) #3

# conclusion: nothing here to pursue

## -------------------------------------------------------------------------
# After have lists, filter on weighted meth difference of exon/intron
## -------------------------------------------------------------------------

head(introns_with_2_cpgs)
intron_data <- output[(output$gene_id %in% introns_with_2_cpgs$gene_id &
                         output$feature == "intron"),]

head(exons_with_2_cpgs)
exon_data <- output[(output$gene_id %in% exons_with_2_cpgs$gene_id &
                         output$feature == "exon"),]

# Read in the weighted methylation levels
weighted_meth <- read_delim("Dcitri_weighted_meth_annotation_by_sex.txt", 
                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
head(weighted_meth)

# Keep only features where the meth level is greatre than background 0.05 in at least one sex
weighted_meth <- weighted_meth[(weighted_meth$female > 0.05 | weighted_meth$male > 0.05),]
weighted_meth <- weighted_meth[!is.na(weighted_meth$chr),]

# remove rows where no data for male or female
weighted_meth <- weighted_meth[!is.na(weighted_meth$male),]
weighted_meth <- weighted_meth[!is.na(weighted_meth$female),]

# Add column for weighted methylation difference between males and females (it's already basically a percentage remember)
weighted_meth$percent_meth_difference_of_feature <- weighted_meth$male - weighted_meth$female

# Keep genes that have at least two diff CpGs in exons or introns and 15%
# feature level difference overall
weighted_meth_exons <- weighted_meth[(weighted_meth$gene_id %in% exon_data$gene_id &
                                       weighted_meth$feature == "exon"),] 
length(unique(weighted_meth_exons$gene_id)) # 39 genes
weighted_meth_exons <- weighted_meth_exons[(weighted_meth_exons$percent_meth_difference_of_feature > 0.15 |
                                             weighted_meth_exons$percent_meth_difference_of_feature < -0.15),]
length(unique(weighted_meth_exons$gene_id)) # 12 genes


weighted_meth_introns <- weighted_meth[(weighted_meth$gene_id %in% intron_data$gene_id &
                                        weighted_meth$feature == "intron"),] 
length(unique(weighted_meth_introns$gene_id)) #50 genes
weighted_meth_introns <- weighted_meth_introns[(weighted_meth_introns$percent_meth_difference_of_feature > 15 |
                                                  weighted_meth_introns$percent_meth_difference_of_feature < -15),] # 31
length(unique(weighted_meth_introns$gene_id)) # 0 genes


# Which sex are these genes hypermethylated in
head(weighted_meth_exons)
weighted_meth_exons$hypermethylated <- "male"
weighted_meth_exons$hypermethylated[weighted_meth_exons$female > weighted_meth_exons$male] <- "female"

female_hyper_exon_genes <- unique(weighted_meth_exons$gene_id[weighted_meth_exons$hypermethylated=="female"]) # 10
length(female_hyper_exon_genes)
male_hyper_exon_genes <-  unique(weighted_meth_exons$gene_id[weighted_meth_exons$hypermethylated=="male"]) # 3
length(male_hyper_exon_genes)
both <- Reduce(intersect, list(female_hyper_exon_genes,male_hyper_exon_genes)) # 1
common_exon <- weighted_meth_exons[weighted_meth_exons$gene_id %in% both,] # two exons in one gene, one female hyper and one male

# Write out a dataframe for later labelling in expression correlation scrips
head(weighted_meth_exons)

write.table(weighted_meth_exons, file="hypermethylated_genes_with_category.txt",
            sep="\t",quote = F, row.names = F, col.names = T)

# Where are the genes in terms of chromosome
weighted_meth_exons$chr[weighted_meth_exons$chr=="DC3.0sc08"] <- "X"
ggplot(weighted_meth_exons, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Differentially Methylated Exons")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))
# Also just by the genes
head(weighted_meth_exons)
gene_data <- weighted_meth_exons[,c(1,3,10)]
gene_data <- gene_data[!duplicated(weighted_meth_exons),]
ggplot(gene_data, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Differentially Methylated Genes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))


## -------------------------------------------------------------------------
# Write out all the gene lists for later use
## -------------------------------------------------------------------------

head(weighted_meth_exons)
write.table(as.data.frame(unique(weighted_meth_exons$gene_id)), file="diff_meth_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)