# Making line and scatter graphs to show relationship of gene exp and meth

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Meth_paired_with_exp")
library(readr)
library(reshape2)
library(Hmisc)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(multcomp)

# -----------------------------------------------
# Read in data
# -----------------------------------------------
# Methylation for genes only
weighted_meth <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation/Dcitri_weighted_meth_genes_only_with_category.txt", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)
diff_meth_genes <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation/hypermethylated_genes_with_category.txt", 
              "\t", escape_double = FALSE, trim_ws = TRUE)

methylation <- merge(weighted_meth, diff_meth_genes, by=c("chr","gene_id"), all=T)
methylation$hypermeth_cat[is.na(methylation$hypermeth_cat)] <- "none"

methylation$male_all <- paste(methylation$male, methylation$male_category, sep="_")
methylation$female_all <- paste(methylation$female, methylation$female_category, sep="_")

methylation <- methylation[,-c(3:6)]
methylation <- reshape2::melt(methylation, id.vars=c("chr","gene_id","hypermeth_cat"))
methylation <- methylation[!duplicated(methylation),]
methylation <- separate(methylation, value, into = c("weighted_meth","meth_level"), sep="_")
colnames(methylation)[4] <- "sex"
methylation$sex <- gsub("_all","", methylation$sex)

# Read in gene expression data
exp_data <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Differential_expression/all_gene_expression_data.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
exp_data <- exp_data[,-c(2:7,11:15)]
exp_data <- reshape2::melt(exp_data, id.vars=c("gene_id","diff_exp","category", "log2FoldChange"))
colnames(exp_data) <- c("gene_id","diff_exp","category","logFC","sex","FPKM")
exp_data$sex <- gsub("_fpkm_mean","", exp_data$sex)

# ----------------------------------------------
# General meth vs general exp full scatter
# -----------------------------------------------
meth_exp <- merge(exp_data, methylation, by=c("gene_id","sex"), all=T) #35934 genes
meth_exp <- meth_exp[!is.na(meth_exp$weighted_meth) & !is.na(meth_exp$FPKM),] #23974 genes
meth_exp$FPKM[meth_exp$FPKM==0] <- 1
meth_exp$logFPKM <- log(meth_exp$FPKM)
meth_exp$logFPKM[meth_exp$logFPKM == -Inf] <- 0

meth_exp$weighted_meth <- as.numeric(meth_exp$weighted_meth)
#write.table(meth_exp, file="all_methylation_expression_data.txt",
 #           col.names = T, row.names = F, quote = F, sep = "\t")

ggplot(meth_exp, aes(x=weighted_meth, y=logFPKM, colour = sex))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


meth_exp_males <- meth_exp[meth_exp$sex=="male",]
ggplot(meth_exp_males, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#44AA99",size=2)+
  geom_smooth(method = "lm",size=2,colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Male")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)

meth_exp_females <- meth_exp[meth_exp$sex=="female",]
ggplot(meth_exp_females, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#DDCC77",size=2)+
  geom_smooth(method = "lm",size=2, colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Female")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)


# Stats
model1<-lm(FPKM~sex*weighted_meth, data=meth_exp)
model2<-lm(FPKM~sex+weighted_meth, data=meth_exp)
anova(model1,model2) # No interaction
summary.lm(model2) # Sex not significiant but meth does predict exp

# ----------------------------------------------
# General meth vs general exp binned scatter
# -----------------------------------------------

meth_exp_females$bins <- as.numeric(cut2(meth_exp_females$weighted_meth , g=100))
meth_exp_males$bins<-as.numeric(cut2(meth_exp_males$weighted_meth , g=100))

female<-as.data.frame(aggregate(meth_exp_females$FPKM, by=list(meth_exp_females$bins), mean))
female$logfpkm <- log10(female$x)
colnames(female)<-c("meth_bin","FPKM","logFPKM")
female$status <- "female"
#plot(female$meth_bin~female$FPKM)

male<-as.data.frame(aggregate(meth_exp_males$FPKM, by=list(meth_exp_males$bins), mean))
male$logfpkm <- log10(male$x)
colnames(male)<-c("meth_bin","FPKM","logFPKM")
male$status <- "male"
#plot(male$meth_bin~male$FPKM)

final_data <- merge(male, female, by="meth_bin")
final_data <- rbind(male, female)

ggplot(data=final_data, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))


# Stats
data1<-reshape2::melt(final_data,id=c("meth_bin","status"))
model1<-lm(value~status*meth_bin, data=data1)
model2<-lm(value~status+meth_bin, data=data1)
anova(model1,model2)# Slight sig interaction
summary.lm(model1) # nothing sig


# ----------------------------------------------
# Binned violin plots (as in Trine's paper)
# -----------------------------------------------
head(meth_exp)

condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

meth_exp$plot_cat <- paste(meth_exp$meth_level, meth_exp$sex, sep="_")
ggplot(meth_exp, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Methylation Category")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("None_female","None_male","Low_female","Low_male","Medium_female",
                            "Medium_male","High_female","High_male"),
                   labels=c("None", "","Low","","Medium","","High",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


# Stats
meth_exp$sex <- as.factor(meth_exp$sex)
meth_exp$meth_level <- as.factor(meth_exp$meth_level)

model1<-lm(FPKM ~ sex * meth_level, data=meth_exp)
model2<-lm(FPKM ~ sex + meth_level, data=meth_exp)
anova(model1,model2) # No sig interaction
summary.lm(model2)

# ----------------------------------------------
# ABOVE but for the X and Autosomes separately
# -----------------------------------------------
meth_exp_X <- meth_exp[meth_exp$chr=="DC3.0sc08",] #1428
meth_exp_A <- meth_exp[!meth_exp$chr=="DC3.0sc08",] #22546

# Scatters
meth_exp_males_X <- meth_exp_X[meth_exp_X$sex=="male",]
meth_exp_males_A <- meth_exp_A[meth_exp_A$sex=="male",]

ggplot(meth_exp_males_A, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#44AA99",size=2)+
  geom_smooth(method = "lm",size=2,colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Male: Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)

meth_exp_females_X <- meth_exp_X[meth_exp_X$sex=="female",]
meth_exp_females_A <- meth_exp_A[meth_exp_A$sex=="female",]

ggplot(meth_exp_females_A, aes(x=weighted_meth, y=logFPKM))+
  geom_point(colour="#DDCC77",size=2)+
  geom_smooth(method = "lm",size=2, colour="black")+
  xlab("Weighted Methylation")+
  ylab("log(FPKM)")+
  ggtitle("Female: Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  xlim(0,1.0)


# Stats
model1<-lm(FPKM~sex*weighted_meth, data=meth_exp_X)
model2<-lm(FPKM~sex+weighted_meth, data=meth_exp_X)
anova(model1,model2) # No interaction
summary.lm(model2) # Nothing signif

# Stats
model1<-lm(FPKM~sex*weighted_meth, data=meth_exp_A)
model2<-lm(FPKM~sex+weighted_meth, data=meth_exp_A)
anova(model1,model2) # No interaction
summary.lm(model2) # Sex not sig but meth exp is


# Binned line
meth_exp_females_X$bins <- as.numeric(cut2(meth_exp_females_X$weighted_meth , g=100))
meth_exp_males_X$bins<-as.numeric(cut2(meth_exp_males_X$weighted_meth , g=100))

meth_exp_females_A$bins <- as.numeric(cut2(meth_exp_females_A$weighted_meth , g=100))
meth_exp_males_A$bins<-as.numeric(cut2(meth_exp_males_A$weighted_meth , g=100))

female_X<-as.data.frame(aggregate(meth_exp_females_X$FPKM, by=list(meth_exp_females_X$bins), mean))
female_X$logfpkm <- log10(female_X$x)
colnames(female_X)<-c("meth_bin","FPKM","logFPKM")
female_X$status <- "female"

female_A<-as.data.frame(aggregate(meth_exp_females_A$FPKM, by=list(meth_exp_females_A$bins), mean))
female_A$logfpkm <- log10(female_A$x)
colnames(female_A)<-c("meth_bin","FPKM","logFPKM")
female_A$status <- "female"

male_X<-as.data.frame(aggregate(meth_exp_males_X$FPKM, by=list(meth_exp_males_X$bins), mean))
male_X$logfpkm <- log10(male_X$x)
colnames(male_X)<-c("meth_bin","FPKM","logFPKM")
male_X$status <- "male"

male_A<-as.data.frame(aggregate(meth_exp_males_A$FPKM, by=list(meth_exp_males_A$bins), mean))
male_A$logfpkm <- log10(male_A$x)
colnames(male_A)<-c("meth_bin","FPKM","logFPKM")
male_A$status <- "male"

final_data_X <- merge(male_X, female_X, by="meth_bin")
final_data_X <- rbind(male_X, female_X)

final_data_A <- merge(male_A, female_A, by="meth_bin")
final_data_A <- rbind(male_A, female_A)

ggplot(data=final_data_A, aes(x=meth_bin, y=logFPKM, colour=status))+
  geom_point(size=2)+
  geom_smooth(method="loess", size=2)+
  xlab("Methylation Rank (Low to High)")+
  ylab("log(FPKM)")+
  ggtitle("Autosomes")+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  scale_colour_manual(breaks = c("female","male"),labels=c("Female","Male"),
                      values=c("#DDCC77","#44AA99"))


# Stats
data1<-reshape2::melt(final_data_X,id=c("meth_bin","status"))
model1<-lm(value~status*meth_bin, data=data1)
model2<-lm(value~status+meth_bin, data=data1)
anova(model1,model2)# no sig interaction
summary.lm(model2) # sig negative relationship...

data1<-reshape2::melt(final_data_A,id=c("meth_bin","status"))
model1<-lm(value~status*meth_bin, data=data1)
model2<-lm(value~status+meth_bin, data=data1)
anova(model1,model2)# no sig interaction
summary.lm(model2) # nothing sig



# Binned violin plots (as in Trine's paper)
head(meth_exp_A)
head(meth_exp_X)

meth_exp_X$plot_cat <- paste(meth_exp_X$meth_level, meth_exp_X$sex, sep="_")
ggplot(meth_exp_A, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Methylation Category")+
  ylab("log(FPKM)")+
  ggtitle("Autosomes")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("None_female","None_male","Low_female","Low_male","Medium_female",
                            "Medium_male","High_female","High_male"),
                   labels=c("None", "","Low","","Medium","","High",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


# Stats
meth_exp_X$sex <- as.factor(meth_exp_X$sex)
meth_exp_X$meth_level <- as.factor(meth_exp_X$meth_level)

model1<-lm(FPKM ~ sex * meth_level, data=meth_exp_X)
model2<-lm(FPKM ~ sex + meth_level, data=meth_exp_X)
anova(model1,model2) # No sig interaction
summary.lm(model2) Nope
