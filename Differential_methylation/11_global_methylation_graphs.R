## -------------------------------------------------------------------------
## Making Fancy Genome-Methylation Differences Graphs
## -------------------------------------------------------------------------

# Load packages etc.
setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation")
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------

# Read in the subsetted file from methylkit 
objectmethbase1 <- read_delim("F_vs_M_objectmethbase.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

objectmethbase1$chrBase <- paste(objectmethbase1$chr, ".", objectmethbase1$start, sep="")
objectmethbase1 <- objectmethbase1[,-3]
objectmethbase1$strand <- gsub(pattern = "\\+", replace = "F", objectmethbase1$strand)


## -------------------------------------------------------------------------

# Setset out each sample
a <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs1","numTs1")]
b <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs2","numTs2")]
c <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs3","numTs3")]
d <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs4","numTs4")]
e <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs5","numTs5")]
f <- objectmethbase1[,c("chrBase","chr","start","strand","coverage1","numCs6","numTs6")]

## -------------------------------------------------------------------------

# Write out each sample
all_files <- list(a,b,c,d,e,f)

for(i in seq_along(all_files)){
  colnames(all_files[[i]])[c(3,5,6,7)] <- c("base","coverage","numCs","numTs")
  all_files[[i]]$freqC <- round((all_files[[i]]$numCs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]]$freqT <- round((all_files[[i]]$numTs/all_files[[i]]$coverage)*100,
                                digits=2)
  all_files[[i]] <- all_files[[i]][-c(6,7)]
  myfile <- file.path("./", paste0(i,"_","subsetted_final.txt"))
  write.table(all_files[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}

# Urgh need to go and rename each file with sample name and condition,
# check the methylkit scritps for the order of samples

## -------------------------------------------------------------------------
# Use methylkit to get the data all togther and plottable 
file.listA <- list("F1_subsetted_final.txt","F2_subsetted_final.txt","F3_subsetted_final.txt",
                   "M1_subsetted_final.txt",
                   "M2_subsetted_final.txt","M3_subsetted_final.txt")

sample_list <- list("F1", "F2", "F3", "M1", "M2", "M3")

# Make a list of all genotypes in the right order = genotype_list

raw_data <- methRead(file.listA,
                     sample.id = sample_list,
                     treatment = c(0,0,0,1,1,1),
                     assembly="Dcitriv3", 
                     context="CpG")

meth_all_data <- unite(raw_data)


## -------------------------------------------------------------------------
# Make a nice dendogram

clusterSamples(meth_all_data, dist="correlation", method="ward", plot=TRUE)

hc <- clusterSamples(meth_all_data, dist="correlation", method="ward", plot=FALSE)

data1 <- dendro_data(hc)
labs <- label(data1)

labs$Sex<- c(rep("Female", 3), rep("Male", 3))

ggplot(segment(data1)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1,
               show.legend = F)+
  geom_text(data=labs,
            aes(label=label, x=x, y=-0.009, colour=Sex,fontface="bold"),
            show.legend = T, angle = 90, size = 5, hjust = 0.5)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text = element_blank(),
        legend.position = "none")+
  scale_colour_manual(values=c("#DDCC77","#44AA99"))

## -------------------------------------------------------------------------
# Make a nice PCA

PCA_data <- PCASamples(meth_all_data, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

PCA_data1$Sex <- c(rep("Female", 3), rep("Male", 3))


percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )


ggplot(PCA_data1, aes(PC1, PC2, colour=Sex))+
  geom_point(size=14)+
  geom_text_repel(aes(label=sample), size=12,show.legend=FALSE, 
                  point.padding = 2, box.padding = 1)+
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.title=element_blank())+
  scale_colour_manual(values=c("#DDCC77","#44AA99"))

# Include a coefficent of variation to measure the dispersion across individual CpGs in males and females
# Function from: https://rcompanion.org/rcompanion/c_02.html
summary.list = function(x)list(
  N.with.NA.removed= length(x[!is.na(x)]),
  Count.of.NA= length(x[is.na(x)]),
  Mean=mean(x, na.rm=TRUE),
  Median=median(x, na.rm=TRUE),
  Max.Min=range(x, na.rm=TRUE),
  # Range=max(Data$ Fish, na.rm=TRUE) - min(Data$ Fish, na.rm=TRUE),
  Variance=var(x, na.rm=TRUE),
  Std.Dev=sd(x, na.rm=TRUE),
  Coeff.Variation.Prcnt=sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100,
  Std.Error=sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)])),
  Quantile=quantile(x, na.rm=TRUE)
)

# Should probably do it by sample
head(objectmethbase1)
objectmethbase1$female_1 <- 1 - ((objectmethbase1$coverage1 - objectmethbase1$numCs1) / objectmethbase1$coverage1)
objectmethbase1$female_2 <- 1 - ((objectmethbase1$coverage2 - objectmethbase1$numCs2) / objectmethbase1$coverage2)
objectmethbase1$female_3 <- 1 - ((objectmethbase1$coverage3 - objectmethbase1$numCs3) / objectmethbase1$coverage3)
objectmethbase1$male_1 <- 1 - ((objectmethbase1$coverage4 - objectmethbase1$numCs4) / objectmethbase1$coverage4)
objectmethbase1$male_2 <- 1 - ((objectmethbase1$coverage5 - objectmethbase1$numCs5) / objectmethbase1$coverage5)
objectmethbase1$male_3 <- 1 - ((objectmethbase1$coverage6 - objectmethbase1$numCs6) / objectmethbase1$coverage6)

summary.list(objectmethbase1$female_1) # 109.0533
summary.list(objectmethbase1$female_2) # 109.0771
summary.list(objectmethbase1$female_3) # 106.7847
summary.list(objectmethbase1$male_1) # 109.9788
summary.list(objectmethbase1$male_2) # 109.4531
summary.list(objectmethbase1$male_3) # 110.5406

## -------------------------------------------------------------------------
# Make a plot to show the methylation level overall for males and females
# Levels taken from the bismark calculation minus the lambda conversion efficienct
# NOTE: it's so similar the graph tells us nothing
library(reshape2)

Female <- c(0.30,0.31,0.31)
Male <- c(0.30,0.30,0.30)

levels <-data.frame(Female, Male)
levels2 <- melt(levels)
levels2 <- levels2[!is.na(levels2$value),]
colnames(levels2) <- c("Sex","Methylation")

ggplot(levels2, aes(y=Methylation, x=Sex,fill=Sex))+
  geom_boxplot()+
  theme_bw()+
  xlab("Sex") +
  ylab("CpG Methylation Level") +
  ylim(0.29,0.32)+
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=30),
        legend.text=element_text(size=30),
        legend.position = "none")+
  scale_fill_manual(values=c("#DDCC77","#44AA99"))
head(levels2)


