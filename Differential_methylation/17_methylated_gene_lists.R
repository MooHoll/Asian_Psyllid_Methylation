## -------------------------------------------------------------------------
# Getting methylated gene lists
## -------------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
## -------------------------------------------------------------------------
annotation <- read_delim("Dcitri_weighted_meth_annotation_by_sex.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
annotation$feature <- as.factor(annotation$feature)

## -------------------------------------------------------------------------
#Genes with meth genes
## -------------------------------------------------------------------------
genes <- annotation[annotation$feature == "gene",]
genes_males <- aggregate(male ~ gene_id + chr, data=genes,FUN = mean)
genes_females <- aggregate(female ~ gene_id + chr, data=genes,FUN = mean)

genes_males$male_category[genes_males$male == 0] <- "None"
genes_males$male_category[genes_males$male > 0 &
                            genes_males$male < 0.3] <- "Low"
genes_males$male_category[genes_males$male >= 0.3 &
                            genes_males$male < 0.7] <- "Medium"
genes_males$male_category[genes_males$male >= 0.7] <- "High"
table(genes_males$male_category)
#high    low medium   none 
#284  15564   1044    703

genes_females$female_category[genes_females$female == 0] <- "None"
genes_females$female_category[genes_females$female > 0 &
                                genes_females$female < 0.3] <- "Low"
genes_females$female_category[genes_females$female >= 0.3 &
                                genes_females$female < 0.7] <- "Medium"
genes_females$female_category[genes_females$female >= 0.7] <- "High"
table(genes_females$female_category)
#high    low medium   none 
#321  15728   1024    609 

head(genes)
# Make a file of general meth genes as a background GO set, meth genes = genes greater than lambda weighted meth 
meth_genes <- genes[genes$female > 0.05 | genes$male > 0.05,]
length(unique(meth_genes$gene_id)) #4429

write.table(as.data.frame(meth_genes$gene_id), file="methylated_genes.txt", sep="\t", quote = F,
            row.names = F, col.names = T)
## -------------------------------------------------------------------------
# Need to pull out the gene list for each and the total gene list as background for the 
# GO enrichment analysis
## -------------------------------------------------------------------------
head(genes_males)
head(genes_females)
all <- merge(genes_males, genes_females, by=c("gene_id","chr"))

write.table(all, file="Dcitri_weighted_meth_genes_only_with_category.txt",
            sep="\t",quote = F, row.names = F, col.names = T)


write.table(unique(genes_males$gene_id[genes_males$male_category=="None"]), 
            file="no_male_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_males$gene_id[genes_males$male_category=="Low"]), 
            file="low_male_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_males$gene_id[genes_males$male_category=="Medium"]), 
            file="medium_male_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_males$gene_id[genes_males$male_category=="High"]), 
            file="high_male_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)

write.table(unique(genes_females$gene_id[genes_females$female_category=="None"]), 
            file="no_female_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_females$gene_id[genes_females$female_category=="Low"]), 
            file="low_female_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_females$gene_id[genes_females$female_category=="Medium"]), 
            file="medium_female_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)
write.table(unique(genes_females$gene_id[genes_females$female_category=="High"]), 
            file="high_female_genes.txt",
            sep="\t",quote = F, row.names = F, col.names = F)

## -------------------------------------------------------------------------
# Figure out if any of these gene categories are enriched on any specific chromosome
## -------------------------------------------------------------------------
genes_males$chr[genes_males$chr == "DC3.0sc08"] <- "X"

ggplot(genes_males, aes(x=chr, fill=factor(male_category, 
                                           levels = c("High", "Medium", "Low","None"))))+
  geom_bar(position="fill")+
  guides()+
  xlab("Chromosome")+
  ylab("Proportion of Genes")+
  ggtitle("Males")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("High", "Medium", "Low", "None"),
                    values=c("#D55E00","#E69F00", "#009E73","#999999"))

genes_females$chr[genes_females$chr == "DC3.0sc08"] <- "X"

ggplot(genes_females, aes(x=chr, fill=factor(female_category, 
                                           levels = c("High", "Medium", "Low","None"))))+
  geom_bar(position="fill")+
  guides()+
  xlab("Chromosome")+
  ylab("Proportion of Genes")+
  ggtitle("Females")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("High", "Medium", "Low", "None"),
                    values=c("#D55E00","#E69F00", "#009E73","#999999"))

