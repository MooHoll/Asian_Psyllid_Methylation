# --------------------------------------------------------------------------
# Identification of chromosome outliers to see if we have X0 or XY system
# --------------------------------------------------------------------------

# Following Kamil's pipeline, refernce: (ask if there is a paper to ref for this pipeline)
# https://github.com/RossLab/PGE/blob/master/sex_chromosome_p_citry.md

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Identification_sex_chroms/windows_for_y")

library(readr)
library(ggplot2)
library(reshape2)

female_depth_file <- read.table('female_.regions.bed', col.names = c('scf', 'start', "end","name","mean_depth"))  
male_depth_file <- read.table('male_.regions.bed', col.names = c('scf', 'start', "end","name","mean_depth"))

female_depth_file$mean_depth[female_depth_file$mean_depth < 1] <- 1
male_depth_file$mean_depth[male_depth_file$mean_depth < 1] <- 1

# Get the average coverage
male_haploid_cov <- mean(male_depth_file$mean_depth) # 49.73
female_haploid_cov <- mean(female_depth_file$mean_depth) # 49.68

# normalization by mean coverage
male_depth_file$male_copies <- male_depth_file$mean_depth / male_haploid_cov
female_depth_file$female_copies <- female_depth_file$mean_depth / female_haploid_cov

all_data <- merge(female_depth_file, male_depth_file, by="name")
all_data <- all_data[,-c(2,3,4)]
colnames(all_data) <- c("name","female_depth","female_copies","scf","start","end","male_depth","male_copies")

all_data$male2female_logratio <- log2(all_data$male_copies / all_data$female_copies)
all_data$scf <- gsub("DC3.0sc", "", all_data$scf)

full_plot <- ggplot(all_data, aes (x=male2female_logratio, fill=scf))+
  geom_histogram(bins=100)+
  xlim(-1.2, 1.2)+
  theme_bw()

full_plot #Doesn't initially appear to show a diverged Y

full_plot +
  xlab("log2(male copies / female copies)")+
  ylab("Count")+
  scale_fill_discrete(name="Chromosome")+
  theme(legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=26))


# Check each chromosome individually
full_plot + facet_grid(scf ~ .) +
  xlab("log2(male copies / female copies)")+
  ylab("Count")+
  theme(legend.position = "none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=22)) # Still looks like no divergeed Y

# Just to be sure with this method I will check each chromosome in it's own plot (more resolution)
for(i in unique(all_data$scf)){
  plot1 <- ggplot(aes(x = male2female_logratio), data = subset(all_data, scf == i)) +
    geom_histogram(bins=100)+
    xlim(-1.2, 1.2)+
    theme_bw()
  ggsave(filename = sprintf('%s.png', i), plot = plot1)
} # Still no diverged Y

all_data_subset <- all_data[all_data$scf=="08",]
ggplot(all_data_subset,aes(x = male2female_logratio)) +
  geom_histogram(bins=100)+
  annotate("text", x=0.5, y=225, label="Chromosome 08", size=10)+
  xlim(-1.2, 1.2)+
  theme_bw()+
  xlab("log2(male copies / female copies)")+
  ylab("Count")+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=26))

# Take a look at the raw sequences
look_for_about0_in_female <- all_data[all_data$female_depth <=3,]

ggplot(look_for_about0_in_female, aes (x=male2female_logratio, fill=scf))+
  geom_histogram(bins=50)+
  xlim(-1.2, 1.2)+
  theme_bw()

depth_data <- look_for_about0_in_female[,c(2,4,7)]
depth_data <- melt(depth_data)

# Nothing with no coverage in females that seems to have coverage in males
ggplot(depth_data, aes (x=value, fill=variable))+
  geom_histogram(bins=50)+
  theme_bw()

