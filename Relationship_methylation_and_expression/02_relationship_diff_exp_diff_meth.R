# -----------------------------------------------
# Relationship of differential methylation and Diff exp
# -----------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Meth_paired_with_exp")

library(readr)
library(ggplot2)
library(FSA)
library(reshape2)
library(Hmisc)
library(tidyr)
library(multcomp)
library(UpSetR)
library(stringr)

# -----------------------------------------------
# Read in all data
# -----------------------------------------------
all_methylation_expression_data <- read_delim("all_methylation_expression_data.txt", 
                                              "\t", escape_double = FALSE, trim_ws = TRUE)

# Add in the weighted methylation difference overall
head(all_methylation_expression_data)
weighted_meth <- all_methylation_expression_data[,c(1,2,9)]
weighted_meth <- weighted_meth[!duplicated(weighted_meth),]
weighted_meth <- spread(weighted_meth, sex, weighted_meth)
weighted_meth$meth_diff <- as.numeric(weighted_meth$male - weighted_meth$female)
weighted_meth <- weighted_meth[,c(1,4)]

all <- merge(all_methylation_expression_data, weighted_meth, by = "gene_id", all=T)

# -----------------------------------------------
# Scatter plots of all differential data
# -----------------------------------------------

# Remove one non-significant outlier messing up the graph
all$logFC <- as.numeric(all$logFC)
all <- all[!(all$logFC < -15),]

all$diff_meth_cat <- "yes"
all$diff_meth_cat[all$hypermethylated=="none"] <- "no"

#Lets take a look, scatter plot all data
ggplot(all, aes(x=meth_diff, y=logFC, colour=diff_exp))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  theme_bw()+
  ggtitle("Differentially Expressed Genes")+
  scale_colour_manual("", breaks=c("no", "yes"),
                      values = c("black","red"),
                      labels= c("No", "Yes"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)

ordered <- all[order(all$diff_meth_cat),]
ggplot(ordered, aes(x=meth_diff, y=logFC, colour=diff_meth_cat))+
  geom_point()+
  xlab("Weighted Methylation Difference")+
  ylab("log(Fold-Change)")+
  theme_bw()+
  ggtitle("Differentially Methylated Genes")+
  scale_colour_manual("", breaks=c("no", "yes"),
                      values = c("black","red"),
                      labels= c("No", "Yes"))+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())+
  geom_vline(xintercept=0.15)+
  geom_vline(xintercept=-0.15)+
  geom_hline(yintercept = 1.5)+
  geom_hline(yintercept = -1.5)+
  xlim(-0.5,0.6)

# -----------------------------------------------
# Binned diff/non-diff methylated gene, violin plots
# -----------------------------------------------

condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Exp levels of diff methylated genes
all$plot_cat <- paste(all$hypermethylated, all$sex, sep="_")
ggplot(all, aes(x=plot_cat, y=logFPKM))+
  geom_violin(aes(fill=sex))+
  xlab("Differential Methylation Category")+
  ylab("log(FPKM)")+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("none_female","none_male","female_female","female_male",
                            "male_female","male_male"),
                   labels=c("None", "None","Female Hypermethylated \nExon","", "Male Hypermethylated \nExon",""))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=16, hjust=0.01),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())


# Stats
all$sex <- as.factor(all$sex)
all$hypermethylated <- as.factor(all$hypermethylated)

model1<-lm(FPKM ~ sex * hypermethylated, data=all)
model2<-lm(FPKM ~ sex + hypermethylated, data=all)
anova(model1,model2) # No sig interaction
summary.lm(model2) # None

# -----------------------------------------------
# Binned diff/non-diff exp gene, violin plots
# -----------------------------------------------

all$log_meth <- log(all$weighted_meth)

all$higher_cat <- "unbiased"
all$higher_cat[all$category=="male_biased"] <- "male_biased"
all$higher_cat[all$category=="male_limited"] <- "male_biased"
all$higher_cat[all$category=="male_biased_extreme"] <- "male_biased"
all$higher_cat[all$category=="female_biased"] <- "female_biased"
all$higher_cat[all$category=="female_limited"] <- "female_biased"
all$plot_column <- paste0(all$sex,"_",all$higher_cat)

ggplot(all, aes(x=plot_column, y=log_meth)) + 
  geom_violin( aes(fill=sex))+
  stat_summary(fun.data=condifence_intervals, color="red",size=1)+
  xlab("Differential Expression Category")+
  ylab("log(Weighted Methylation)")+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("female_unbiased","male_unbiased", "female_female_biased",
                            "male_female_biased","female_male_biased","male_male_biased"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Stats
model1<-lm(weighted_meth ~ sex * higher_cat, data=all)
model2<-lm(weighted_meth ~ sex + higher_cat, data=all)
anova(model1,model2) # No sig interaction
summary.lm(model2) # Sig: unbiased genes have higher methylation than biased genes

# -----------------------------------------------
# Same binned diff exp but as bar plots
# -----------------------------------------------

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
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summary_data<-summarySE(all, measurevar = "weighted_meth", 
                        groupvars = c("plot_column","sex"))

ggplot(summary_data, aes(x=plot_column, y=weighted_meth, fill=sex))+
  geom_bar(position = position_dodge(), stat = "identity")+
  geom_errorbar(aes(ymin=weighted_meth-ci, ymax=weighted_meth+ci),
                width=.2,
                position = position_dodge(.9))+
  xlab("Differentially Expressed Category")+
  ylab("log(Average Weighted Methylation)")+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("female_unbiased","male_unbiased", "female_female_biased",
                            "male_female_biased","female_male_biased","male_male_biased"),
                   labels=c("Unbiased", "Unbiased",
                            "Female\nBiased", "Female\nBiased",
                            "Male\nBiased","Male\nBiased"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Confirm there are no genes which are both diff exp and diff meth
look <- all[!(all$diff_exp=="no") & !(all$hypermeth_cat=="none"),] # 0


