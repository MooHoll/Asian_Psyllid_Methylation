## -------------------------------------------------------------------------
# Re-make meth over feature graph for diff levels of methylation
## -------------------------------------------------------------------------

# NOTE: We want to split exons into groups if 1-3 are say more methylated
# come back to this script once we have figured this out

# Also can plot TE methylation by types of TE and can plot meth levels over chromosomes


setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation")
library(readr)
library(doBy)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(scales)
library(ggpubr)
## -------------------------------------------------------------------------
annotation <- read_delim("Dcitri_weighted_meth_annotation_by_sex.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

annotation <- annotation[,-c(4:6)]

melted_annot <- melt(annotation, id.vars = c("chr","feature","gene_id") )
colnames(melted_annot) <- c("Chromosome","Feature","ID","Sex","Weighted_Methylation")

# Remove rows where NA in one sex
melted_annot <- melted_annot[!is.na(melted_annot$Weighted_Methylation),]

# Only a couple of splice sites so remove as well
melted_annot <- melted_annot[!melted_annot$Feature == "non_canonical_five_prime_splice_site",]
melted_annot <- melted_annot[!melted_annot$Feature == "non_canonical_three_prime_splice_site",]

# Remove mRNA and CDS not really informative
melted_annot <- melted_annot[!melted_annot$Feature == "CDS",]
melted_annot <- melted_annot[!melted_annot$Feature == "mRNA",]

## -------------------------------------------------------------------------
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
## -------------------------------------------------------------------------
# Normal graph with all information

summary_all<-summarySE(melted_annot, measurevar = "Weighted_Methylation", 
                       groupvars = c("Feature","Sex"))

ggplot(summary_all, aes(x=Feature, y=Weighted_Methylation, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

# Stats
library(multcomp)
head(melted_annot)
model1<-lm(Weighted_Methylation ~ Sex * Feature , data=melted_annot)
model2<-lm(Weighted_Methylation ~ Sex + Feature , data=melted_annot)
anova(model1,model2) # No interaction
summary.lm(model2) # Everything sig

for_stats <- melted_annot
for_stats$SHD<-interaction(for_stats$Sex,
                           for_stats$Feature)
model1_new<-lm(Weighted_Methylation~-1+SHD, data=for_stats)
summary(glht(model1_new,linfct=mcp(SHD="Tukey"))) 

## -------------------------------------------------------------------------
# Different categories
head(melted_annot)
#limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))
look <- melted_annot[melted_annot$Feature=="intergenic",]
hist(look$Weighted_Methylation[look$Weighted_Methylation > 0.05], main = "Intergenic methylation level distribution",
     xlab = "Weighted Methylation", cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)


## -------------------------------------------------------------------------
# Frequency of features being high/low methylated
head(melted_annot)
meth_low <- melted_annot[melted_annot$Weighted_Methylation < 0.3,]
meth_medium <- melted_annot[melted_annot$Weighted_Methylation > 0.3 &
                              melted_annot$Weighted_Methylation < 0.7,]
meth_high <- melted_annot[melted_annot$Weighted_Methylation > 0.7,]
meth_none <- melted_annot[melted_annot$Weighted_Methylation ==0,]

melted_annot$bins<-"low"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.3 &
                    melted_annot$Weighted_Methylation < 0.7] <-"medium"
melted_annot$bins[melted_annot$Weighted_Methylation > 0.7] <-"high"
melted_annot$bins[melted_annot$Weighted_Methylation ==0] <-"none"

melted_annot$combined <- paste0(melted_annot$Feature, "_", melted_annot$Sex)
melted_meth_stuff_2 <- melted_annot[,-c(3,5)]
melted_meth_stuff_2$counts <- with(melted_meth_stuff_2, 
                                  ave(bins, Feature, Sex, bins, FUN=length))
plot_data <- melted_meth_stuff_2[!duplicated(melted_meth_stuff_2),]

plot_data_prom <- subset(plot_data, bins =="low")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b1<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("Genomic Feature")+
#  ylab("Count")+
  ggtitle("Low Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="medium")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b2<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  #xlab("Genomic Feature")+
  #ylab("Count")+
  ggtitle("Medium Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
 # xlab("")+
  ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title.y=element_text(size=12),
        axis.title.x = element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  xlab("Genomic Feature")+
  ylab("Count")+
 # ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

## -------------------------------------------------------------------------
plot_data_prom <- subset(plot_data, bins =="none")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b4 <- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("Genomic Feature")+
  #  ylab("Count")+
  ggtitle("No Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))

plot_data_prom <- subset(plot_data, bins =="high")
plot_data_prom$counts <- as.numeric(plot_data_prom$counts)

b3_1<- ggplot(plot_data_prom, aes(x=Feature, fill=Sex, y=counts))+
  geom_bar(position = position_dodge(), stat = "identity")+
  theme_bw()+
  # xlab("")+
 # ylab("Count")+
  ggtitle("High Methylation")+
  theme(axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=10),
        axis.title=element_blank(),
        plot.title=element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12))+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(breaks = c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"),
                   labels = c("Promoter","5' UTR","3' UTR","Gene","Exon","Intron","TE","Intergenic"),
                   limits =c("promoter","five_prime_UTR","three_prime_UTR","gene","exon","intron","TE","intergenic"))


## -------------------------------------------------------------------------
# Slightly diff figures
levels_across_feature <- ggarrange(b3_1,b2,b1,b4,
                                   ncol=2, nrow=2, common.legend = TRUE, legend="right")

annotate_figure(levels_across_feature, 
                 left = text_grob("Count", 
                                 color = "black", rot = 90, size=14),
                bottom = text_grob("Genomic Feature", 
                                   color = "black", size =12,
                                   hjust = 0.75 ))


## -------------------------------------------------------------------------
# Take a look at different TE cetegories
head(melted_annot)

tes <- melted_annot[melted_annot$Feature=="TE",]
head(tes)
unique(tes$ID)

summary_all<-summarySE(tes, measurevar = "Weighted_Methylation", 
                       groupvars = c("ID","Sex"))

ggplot(summary_all, aes(x=ID, y=Weighted_Methylation, fill=Sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=Weighted_Methylation-ci, ymax=Weighted_Methylation+ci),
                width=.2,
                position = position_dodge(.9))+
  theme_bw()+
  xlab("Transposable Element Type")+
  ylab("Weighted Methylation Level")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))

# Let's also just count up the TE content generally
annotation_rawdata <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files/Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters_plusTEs_plusIntergenic.txt", 
                                                                                             "\t", escape_double = FALSE, trim_ws = TRUE)
tes <- annotation_rawdata[annotation_rawdata$feature == "TE",]
head(tes)
tes$chr_cat <- "A"
tes$chr_cat[tes$chr=="DC3.0sc08"] <- "X"

ggplot(tes, aes(x=gene_id, fill=chr_cat))+
  geom_bar(position = position_dodge2(preserve = "single", padding = 0), stat = "count")+
  theme_bw()+
  xlab("Transposable Element Type")+
  ylab("Count")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

tes$size <- tes$end - tes$start
head(tes)

# Need total genome size - X chromosome size here
genome_fai <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files/Diaci_v3.0.ref_nospaces.fa.fai", 
                                         "\t", escape_double = FALSE, col_names = FALSE, 
                                         trim_ws = TRUE)
genome_fai <- genome_fai[,c(1,2)]
# 08 = 28791834
genome_fai <- genome_fai[!genome_fai$X1 =="DC3.0sc08",]
sum(genome_fai$X2) # 445174712

tes$chr_size <- 445174712
tes$chr_size[tes$chr=="DC3.0sc08"] <- 28791834
head(tes)

tes <- tes[,-c(1,2,3,4)]
te_summary <- summaryBy(size ~ chr_cat + gene_id + chr_size,
                        data = tes, FUN=sum)
head(te_summary)
te_summary$percent <- 100*(te_summary$size.sum / te_summary$chr_size)

ggplot(te_summary, aes(x=gene_id, fill=chr_cat, y = percent))+
  geom_bar(position = position_dodge2(preserve = "single", padding = 0), stat = "identity")+
  theme_bw()+
  xlab("Transposable Element Type")+
  ylab("Percentage of Chromosome Occupied")+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(angle=45,hjust=1,size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Also quick look at the total percentage of the genome occupied by TEs
445174712 + 28791834 # 473966546
head(te_summary)

sum(te_summary$size.sum) # 15658460
100*(15658460/473966546) # 3.3%

