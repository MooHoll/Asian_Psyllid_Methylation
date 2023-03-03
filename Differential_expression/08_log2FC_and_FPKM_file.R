#------------------------------------------------------------------
# Making a file with FPKM and log2FC for later incorporation with methylation data
#------------------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Differential_expression")

library(readr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggforce)

#------------------------------------------------------------------
# log2FC from diff exp
diff_exp_log2FC_all <- read_delim("diff_exp_output_log2FC_all_genes.txt", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
diff_exp_log2FC_all <- diff_exp_log2FC_all[,c(2,5,6,7)]
colnames(diff_exp_log2FC_all)[4] <- "gene_id"

#------------------------------------------------------------------
# all samples for FPKM
F1 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_F2-1_FRAS202200995-1r.genes.results", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
F1 <- F1[,c(1,7)]
colnames(F1) <- c("gene_id", "F1_FPKM")


F2 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_F2-2_FRAS202200996-1r.genes.results", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
F2 <- F2[,c(1,7)]
colnames(F2) <- c("gene_id", "F2_FPKM")


F3 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_F2-3_FRAS202200997-1r.genes.results", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
F3 <- F3[,c(1,7)]
colnames(F3) <- c("gene_id", "F3_FPKM")


M1 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_M2-1_FRAS202200998-1r.genes.results", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
M1 <- M1[,c(1,7)]
colnames(M1) <- c("gene_id", "M1_FPKM")


M2 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_M2-2_FRAS202200999-1r.genes.results", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
M2 <- M2[,c(1,7)]
colnames(M2) <- c("gene_id", "M2_FPKM")


M3 <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Counts/trim_M2-3_FRAS202201000-1b.genes.results", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)
M3 <- M3[,c(1,7)]
colnames(M3) <- c("gene_id", "M3_FPKM")

#------------------------------------------------------------------
# Put it all together
fpkm_data <- Reduce(merge, list(F1,F2,F3,M1,M2,M3))

fpkm_data$female_fpkm_mean <- (fpkm_data$F1_FPKM + fpkm_data$F1_FPKM + fpkm_data$F1_FPKM)/3
fpkm_data$male_fpkm_mean <- (fpkm_data$M1_FPKM + fpkm_data$M1_FPKM + fpkm_data$M1_FPKM)/3

all_data <- merge(fpkm_data, diff_exp_log2FC_all, by = "gene_id")

# 12420/19049 annotated genes as removed genes with low count during diff exp pipeline (mean fpkm < 10)

# Measure of specificity for gene expression
all_data$female_sq <- all_data$female_fpkm_mean*all_data$female_fpkm_mean
all_data$male_sq <- all_data$male_fpkm_mean*all_data$male_fpkm_mean
all_data$SPM <- all_data$female_sq/(all_data$female_sq +all_data$male_sq )

ggplot(all_data, aes ( x= SPM))+
  geom_histogram(colour="black", bins=50)+
  xlab("SPM Relative to Females")+
  ylab("Number of Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

#------------------------------------------------------------------
# Add gene categories for plotting
all_data$diff_exp <- "no"
all_data$diff_exp[all_data$padj < 0.05] <- "yes"
all_data$diff_exp[all_data$log2FoldChange < 1.5 & all_data$log2FoldChange > -1.5] <- "no"

table(all_data$diff_exp) # 1259 = yes

all_data$category <- "unbiased"
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange > 1.5] <- "male_biased"
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange < -1.5] <- "female_biased"
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange > 10] <- "male_biased_extreme"
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange < -10] <- "female_biased_extreme" # NONE

table(all_data$category)

all_data$category[all_data$diff_exp =="yes" & all_data$female_fpkm_mean==0] <- "male_limited"
all_data$category[all_data$diff_exp =="yes" & all_data$male_fpkm_mean==0] <- "female_limited"

table(all_data$category)

# Take an eyeball of the male limited genes
head(all_data)
look <- all_data[all_data$category=="male_limited",]
boxplot(look$male_fpkm_mean)
mean(look$male_fpkm_mean)
median(look$male_fpkm_mean)
range(look$male_fpkm_mean)
nrow(look[look$male_fpkm_mean<10,]) #88/484 - ok I'm convinced this is real

# Goodness of fit
observed = c(1164, 95)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 



write.table(all_data, file="all_gene_expression_data.txt", sep="\t", quote=F, col.names = T, row.names = F )
#------------------------------------------------------------------
# Plots
all_data_plot <- all_data[!all_data$category == "unbiased",]
all_data_plot$biased[all_data_plot$diff_exp =="yes" & all_data_plot$log2FoldChange > 1.5] <- "male_biased"
all_data_plot$biased[all_data_plot$diff_exp =="yes" & all_data_plot$log2FoldChange < -1.5] <- "female_biased"
all_data_plot$category <- gsub(".*_","",all_data_plot$category)

ggplot(all_data_plot, aes(x=biased, fill=category))+
  geom_bar()+
  theme_bw()+
  xlab("Gene Expression Category")+
  ylab("Number of Genes")+
  scale_fill_manual("",breaks=c("biased", "extreme","limited"),
                    labels = c("Fold-change > 1.5","Fold-change > 10","Sex-Limited"),
                    values = c("grey75", "grey45","grey14"))+
  theme_bw()+
  scale_x_discrete(limits=c("female_biased", "male_biased"),
                   labels=c("Female Biased", "Male Biased"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))
#------------------------------------------------------------------
head(all_data)
all_data$female_log_fpkm <- log10(all_data$female_fpkm_mean)
all_data$male_log_fpkm <- log10(all_data$male_fpkm_mean)

scatter_data <- all_data[,c(1,17:19)]

ggplot(scatter_data, aes(x=female_log_fpkm, y=male_log_fpkm, colour=category))+
  geom_point(size=3)+
  xlim(-1.5,4)+
  ylim(-1.5,4)+
  scale_colour_manual("",breaks=c("unbiased","male_biased","male_biased_extreme","male_limited",
                               "female_biased","female_limited"),
                      labels=c("Unbiased","Male Biased","Male Extreme Biased","Male Limited",
                               "Female Biased","Female Limited"),
                      values=c("grey","#44AA99","springgreen4","darkgreen","#DDCC77","goldenrod"))+
  theme_bw()+
  xlab("Log10(Female FPKM)")+
  ylab("Log10(Male FPKM)")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

#------------------------------------------------------------------
# Look for X chromosome enrichment of sex biased genes

head(all_data)
all_data1 <- separate(all_data, gene_id, into=c("chromosome", "gene_id"), sep="g")
head(all_data1)

autosomes <- all_data1[!all_data1$chromosome == "Dcitr08",]
x_chorm <- all_data1[(all_data1$chromosome == "Dcitr08"),]

ggplot(autosomes, aes(x=category))+
  geom_bar(stat = "count")

ggplot(x_chorm, aes(x=category))+
  geom_bar(stat = "count")


autosomes <- autosomes[,c(1,18)]
x_chorm <- x_chorm[,c(1,18)]

autosomes_percent <- autosomes %>% 
  group_by(category) %>% 
  dplyr::summarise(count = n()) %>% 
  mutate(perc = 100*(count/sum(count)))

x_chorm_percent <- x_chorm %>% 
  group_by(category) %>% 
  dplyr::summarise(count = n()) %>% 
  mutate(perc = 100*(count/sum(count)))

ggplot(autosomes_percent, aes(x=category, y=perc, fill=category))+
  geom_bar(stat="identity")+
  facet_zoom(ylim = c(0, 8))+
  xlab("Gene Expression Category")+
  ylab("Percentage of Genes (%)")+
  ggtitle("Autosomes")+
  scale_fill_manual("",limits=c("unbiased","male_biased","male_biased_extreme","male_limited",
                                  "female_biased","female_limited"),
                      values=c("darkgrey","#44AA99","springgreen4","darkgreen","#DDCC77","goldenrod"))+
  scale_x_discrete(limits=c("unbiased","male_biased","male_biased_extreme","male_limited",
                            "female_biased","female_limited"),
                   labels=c("Unbiased","Male Biased","Male Extreme Biased","Male Limited",
                            "Female Biased","Female Limited"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))
  

ggplot(x_chorm_percent, aes(x=category, y=perc, fill=category))+
  geom_bar(stat="identity")+
  facet_zoom(ylim = c(0, 8))+
  xlab("Gene Expression Category")+
  ylab("Percentage of Genes (%)")+
  ggtitle("X Chromosome")+
  scale_fill_manual("",limits=c("unbiased","male_biased","male_biased_extreme","male_limited",
                                "female_biased","female_limited"),
                    values=c("darkgrey","#44AA99","springgreen4","darkgreen","#DDCC77","goldenrod"))+
  scale_x_discrete(limits=c("unbiased","male_biased","male_biased_extreme","male_limited",
                            "female_biased","female_limited"),
                   labels=c("Unbiased","Male Biased","Male Extreme Biased","Male Limited",
                            "Female Biased","Female Limited"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=20, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))

x_chorm_percent
autosomes_percent

# Take the percentage sum of female_biased and female_limited (or for male)
# are compare between X and A with an expected split of 50/50

# Female X = 0.945
# Female A = 0.7532

# Goodness of fit
observed = c(0.945, 0.7532)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) # X-squared = 0.021662, df = 1, p-value = 0.883

# Male X = 0.27
# Male A = 9.9525

# Goodness of fit
observed = c(0.27, 9.9525)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) # X-squared = 9.171, df = 1, p-value = 0.002459

