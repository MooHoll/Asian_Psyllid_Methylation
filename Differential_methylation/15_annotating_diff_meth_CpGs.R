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


# Here we are saying a min of 2 diff meth cpgs per feature to count
intron_data <- output[output$feature=="intron",]
number_diff_cpgs_per_intron <- count(intron_data, gene_id)
mean(number_diff_cpgs_per_intron$n) #1-7, mean = 1.349, median = 1
hist(number_diff_cpgs_per_intron$n)
nrow(number_diff_cpgs_per_intron[number_diff_cpgs_per_intron$n > 2,]) #27
introns_with_2_cpgs <- subset(number_diff_cpgs_per_intron, n >2)

exon_data <- output[output$feature=="exon",]
number_diff_cpgs_per_exon <- count(exon_data, gene_id)
median(number_diff_cpgs_per_exon$n) #1-7, mean = 1.256, median = 1
hist(number_diff_cpgs_per_exon$n)
nrow(number_diff_cpgs_per_exon[number_diff_cpgs_per_exon$n > 2,]) #14
exons_with_2_cpgs <- subset(number_diff_cpgs_per_exon, n >2)

UTR3_data <- output[output$feature=="three_prime_UTR",]
number_diff_cpgs_per_UTR3 <- count(UTR3_data, gene_id)
range(number_diff_cpgs_per_UTR3$n) #1-6, mean = 1.5, median = 1
hist(number_diff_cpgs_per_UTR3$n)
nrow(number_diff_cpgs_per_UTR3[number_diff_cpgs_per_UTR3$n > 2,]) #3 
UTR3_with_2_cpgs <- subset(number_diff_cpgs_per_UTR3, n >2)

UTR5_data <- output[output$feature=="five_prime_UTR",]
number_diff_cpgs_per_UTR5 <- count(UTR5_data, gene_id)
median(number_diff_cpgs_per_UTR5$n) #1-2, mean = 1.13, median = 1
hist(number_diff_cpgs_per_UTR5$n)
nrow(number_diff_cpgs_per_UTR5[number_diff_cpgs_per_UTR5$n > 2,]) #0 
UTR5_with_2_cpgs <- subset(number_diff_cpgs_per_UTR5, n >2)

promoter_data <- output[output$feature=="promoter",]
number_diff_cpgs_per_prom <- count(promoter_data, gene_id)
range(number_diff_cpgs_per_prom$n) #1-2, mean = 1.21, median = 1
hist(number_diff_cpgs_per_prom$n)
nrow(number_diff_cpgs_per_prom[number_diff_cpgs_per_prom$n > 2,]) #0 
proms_with_2_cpgs <- subset(number_diff_cpgs_per_prom, n >2)

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

# Filter by 3cpgs per TE
unique_TEs_only <- unique(TEs_subset$uniquie_identified)
number_diff_cpgs_per_TE<- count(TEs_subset,uniquie_identified )
TEs_with_2_cpgs <- subset(number_diff_cpgs_per_TE, n >2) #1

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

# Add column for weighted methylation difference between males and females as a %
weighted_meth$percent_meth_difference_of_feature <- ((weighted_meth$male -
                                                        weighted_meth$female) / weighted_meth$male)*100

# Keep genes that have at least two diff CpGs in exons or introns and 15%
# feature level difference overall
weighted_meth_exons <- weighted_meth[(weighted_meth$gene_id %in% exon_data$gene_id &
                                       weighted_meth$feature == "exon"),] #34 exons
length(unique(weighted_meth_exons$gene_id)) # 13 genes
weighted_meth_exons <- weighted_meth_exons[(weighted_meth_exons$percent_meth_difference_of_feature > 15 |
                                             weighted_meth_exons$percent_meth_difference_of_feature < -15),] # 21
length(unique(weighted_meth_exons$gene_id)) # 11 genes


weighted_meth_introns <- weighted_meth[(weighted_meth$gene_id %in% intron_data$gene_id &
                                        weighted_meth$feature == "intron"),] #50 introns
length(unique(weighted_meth_introns$gene_id)) # 23 genes
weighted_meth_introns <- weighted_meth_introns[(weighted_meth_introns$percent_meth_difference_of_feature > 15 |
                                                  weighted_meth_introns$percent_meth_difference_of_feature < -15),] # 31
length(unique(weighted_meth_introns$gene_id)) # 22 genes

# How many genes in common?
exon_gene_ids <- as.data.frame(unique(weighted_meth_exons$gene_id))
colnames(exon_gene_ids) <- "gene_id"

intron_gene_ids <- as.data.frame(unique(weighted_meth_introns$gene_id))
colnames(intron_gene_ids) <- "gene_id"

both <- Reduce(intersect, list(exon_gene_ids,intron_gene_ids)) # 1 - Dcitr02g05900.1

# Which sex are these genes hypermethylated in
head(weighted_meth_exons)
weighted_meth_exons$hypermethylated <- "male"
weighted_meth_exons$hypermethylated[weighted_meth_exons$female > weighted_meth_exons$male] <- "female"

female_hyper_exon_genes <- unique(weighted_meth_exons$gene_id[weighted_meth_exons$hypermethylated=="female"]) # 7
length(female_hyper_exon_genes)
male_hyper_exon_genes <-  unique(weighted_meth_exons$gene_id[weighted_meth_exons$hypermethylated=="male"]) # 7
length(male_hyper_exon_genes)
both <- Reduce(intersect, list(female_hyper_exon_genes,male_hyper_exon_genes)) # 3
common_exon <- weighted_meth_exons[weighted_meth_exons$gene_id %in% both,] # all have 3 diff exons, all two female hyper and one male hyper

head(weighted_meth_introns)
weighted_meth_introns$hypermethylated <- "male"
weighted_meth_introns$hypermethylated[weighted_meth_introns$female > weighted_meth_introns$male] <- "female"

female_hyper_intron_genes <- unique(weighted_meth_introns$gene_id[weighted_meth_introns$hypermethylated=="female"]) # 14
length(female_hyper_intron_genes)
male_hyper_intron_genes <-  unique(weighted_meth_introns$gene_id[weighted_meth_introns$hypermethylated=="male"]) # 10
length(male_hyper_intron_genes)
both <- Reduce(intersect, list(female_hyper_intron_genes,male_hyper_intron_genes)) # 2
length(both)
common_intron <- weighted_meth_introns[weighted_meth_introns$gene_id %in% both,] 


# Where are the genes in terms of chromosome
ggplot(exon_data, aes(x=chr, fill=hypermethylated))+
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

ggplot(intron_data, aes(x=chr, fill=hypermethylated))+
  geom_bar()+
  guides()+
  xlab("Chromosome")+
  ylab("Number of Differentially Methylated Introns")+
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

head(weighted_meth_introns) 
write.table(as.data.frame(unique(weighted_meth_introns$gene_id)), file="diff_meth_intron_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)



# HERE XXXXXX


# Also forthe below output the weighted meth scores for correlations
both # genes with both exon and promotor 1522
write.table(as.data.frame(both), file="./diff_meth_gene_lists/diff_meth_common_promotor_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

prom_only <- promotor_genes_only[!promotor_genes_only %in% both] # 1187
write.table(as.data.frame(prom_only), file="./diff_meth_gene_lists/diff_meth_unique_promotor_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)

exon_only <- exon_genes_only[!exon_genes_only %in% both] # 1214
write.table(as.data.frame(exon_only), file="./diff_meth_gene_lists/diff_meth_unique_exon_geneIDs.txt", 
            sep="\t", quote = F, col.names = F, row.names = F)


both_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% both,]
write.table(as.data.frame(both_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_common_promotor_exon_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

prom_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% prom_only,]
write.table(as.data.frame(prom_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_unique_promotor_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

exon_with_weighted_meth <- exon_genes[exon_genes$gene_id %in% exon_only,]
write.table(as.data.frame(exon_with_weighted_meth), file="./diff_meth_gene_lists/diff_meth_unique_exon_all_info.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

# Make file of all proms and if they're diff or not for Peter
prom_with_weighted_meth <- promotor_genes[promotor_genes$gene_id %in% promotor_genes_only,]
head(prom_with_weighted_meth)
for_peter <- prom_with_weighted_meth[,c(4,5,6,7,9,10,11)]
head(for_peter2)

for_peter2 <- for_peter[!duplicated(for_peter),]
for_peter2$sig_diff_methylated <- "yes"

head(annotation)
annotation_proms <- subset(annotation, feature=="promotors_2000bp")
annotation_proms <- annotation_proms[,c(1,2,3,4,5,7,8)]
nrow(annotation_proms)

not_sig_proms <- annotation_proms[!annotation_proms$gene_id %in% for_peter2$gene_id,]
head(not_sig_proms)
not_sig_proms$sig_diff_methylated <- "no"

all_data_proms <- rbind(for_peter2, not_sig_proms)
head(all_data_proms)
nrow(all_data_proms)
write.table(as.data.frame(all_data_proms), file="pcitri_diff_meth_promotors.txt", 
            sep="\t", quote = F, col.names = T, row.names = F)

FPKM_values_by_sex <- read_delim("~/Dropbox/Edinburgh/Sex-specific-mealybugs/transcription/FPKM_values_by_sex.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
data_wide <- dcast(FPKM_values_by_sex, gene_id ~ origin, value.var="FPKM")

head(all_data_proms)
all_data_proms$hyper_status <- ifelse(all_data_proms$male_mean_weightedMeth >
                                        all_data_proms$female_mean_weightedMeth,
                                  "male", "female")

all_data_proms$hyper_status[all_data_proms$sig_diff_methylated == "no"] <- "not_sig"

both_data <- merge(all_data_proms, data_wide, by="gene_id")
head(both_data)
both_data <- both_data[,c(9,10,11)]
both_data_melt <- melt(both_data)
head(both_data_melt)

#### Define summary function (ref:http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 

summary_all<-summarySE(both_data_melt, measurevar = "value", 
                       groupvars = c("hyper_status","variable"))
head(summary_all)

ggplot(summary_all, aes(x=hyper_status, y=value, fill=variable))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Gene Set")+
  ylab("FPKM")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("female","male","not_sig"),
                   labels = c("Female Hypermethylated","Male Hypermethylated",
                              "Not Significant"))

both_data <- merge(all_data_proms, data_wide, by="gene_id")
head(both_data)
both_data$hyper_status[(both_data$male_mean_weightedMeth < 0.01 &
                                both_data$female_mean_weightedMeth < 0.01)]<-"not_meth"
both_data <- both_data[,c(9,10,11)]
both_data_melt <- melt(both_data)
head(both_data_melt)

summary_all<-summarySE(both_data_melt, measurevar = "value", 
                       groupvars = c("hyper_status","variable"))
head(summary_all)

ggplot(summary_all, aes(x=hyper_status, y=value, fill=variable))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  ggtitle("Promotors")+
  xlab("Gene Set")+
  ylab("FPKM")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("pink1","steelblue1"))+
  scale_x_discrete(breaks = c("female","male","not_sig","not_meth"),
                   labels = c("Female\nHypermethylated","Male\nHypermethylated",
                              "Not Significant","Unmethylated"),
                   limits = c("not_meth","not_sig","female","male"))
