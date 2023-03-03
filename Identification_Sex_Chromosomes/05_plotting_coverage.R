# --------------------------------------------------------------------------
# Identification of chromosome outliers to see if we have X0 or XY system
# --------------------------------------------------------------------------

# Following Kamil's pipeline, refernce: (ask if there is a paper to ref for this pipeline)
# https://github.com/RossLab/PGE/blob/master/sex_chromosome_p_citry.md

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Identification_sex_chroms")

library(readr)
library(ggplot2)
library(reshape2)

female_depth_file <- read.table('female_depth.txt', col.names = c('scf', 'female_depth'))  
male_depth_file <- read.table('male_depth.txt', col.names = c('scf', 'male_depth'))

# Had an issue getting the chromosome lengths as samtools faidx would not
# work because the sequences containsed N's on different lines, urgh
# Found this: 
# awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' Diaci_v3.0.ref.fa | uniq > Dcitri_chroms_length.txt
# here: https://www.biostars.org/p/373962/
# which works! Happy times!

reference_index <-read_table2("Dcitri_chroms_length.txt", col_names = FALSE)
colnames(reference_index) <- c("scf","len")
# remove the unplaced scaffolds from the ref index file
reference_index <- reference_index[-14,]

# preparing and merging of the three data datables
depth_tab <- merge(merge(reference_index, female_depth_file), male_depth_file)

# calculate averge coverage (depth normalized by scf length -> maybe should be redone with non-0 covered nts)
depth_tab$female_depth <- depth_tab$female_depth / depth_tab$len
depth_tab$male_depth <- depth_tab$male_depth / depth_tab$len

male_haploid_cov <- mean(depth_tab$male_depth) # 50.5
female_haploid_cov <- mean(depth_tab$female_depth) # 51.2 

# normalization by mean coverage
male_copies <- depth_tab$male_depth / male_haploid_cov
female_copies <- depth_tab$female_depth / female_haploid_cov

male2female_logratio <- as.data.frame(log2(male_copies / female_copies))
colnames(male2female_logratio) <- "ratio"

ggplot(male2female_logratio, aes (x=ratio))+
  geom_histogram()+
  xlim(-0.8, 0.8)

# plot for the coverage / mean coverage
depth_tab$male_copies <-depth_tab$male_depth/male_haploid_cov
depth_tab$female_copies <-depth_tab$female_depth/female_haploid_cov

plot_depth_data <- depth_tab[,c(1,5,6)]
plot_depth_data2 <- melt(plot_depth_data)
colnames(plot_depth_data2) <- c("Chromosome","Sex","Copies")

plot_depth_data2$Chromosome <- gsub("DC3.0sc", "", plot_depth_data2$Chromosome)

ggplot(plot_depth_data2, aes(x=Chromosome, y= Copies, fill=Sex))+
  geom_bar(stat="identity", position="dodge", colour="black")+
  scale_fill_manual(values = c("#44AA99","#DDCC77"),
                    limits = c("male_copies","female_copies"),
                    labels = c("Male","Female"))+
  theme_bw()+
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=20),
        legend.text = element_text(size=28),
        legend.title = element_blank())

# 6699CC
#   scale_fill_manual(values = c("#0072B2","#E69F00"),
