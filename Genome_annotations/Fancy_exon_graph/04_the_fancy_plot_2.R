## -------------------------------------------------------------------------
# Filtering diff meth CpGs from methylkit with weighted meth of features
## -------------------------------------------------------------------------
setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Genome_wide")

library(readr)
library(doBy)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)

female_weighted_meth <- read_delim("female_weighted_meth.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
male_weighted_meth <- read_delim("male_weighted_meth.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)

head(female_weighted_meth)
unique(female_weighted_meth$id)

# Cut down the dataframes for easier use
female_weighted_meth <- female_weighted_meth[,c(2,8)]
female_weighted_meth$sex <- "female"
male_weighted_meth <- male_weighted_meth[,c(2,8)]
male_weighted_meth$sex <- "male"

# Take the average methylation level per bin across all features
female_summary <- summaryBy(weightedMeth ~ id + sex, data = female_weighted_meth, FUN=mean)
male_summary <- summaryBy(weightedMeth ~ id + sex, data = male_weighted_meth, FUN=mean)

all_data <- rbind(female_summary, male_summary, by = "id")

# Add new column to show the bins for easier plotting (maybe?)
all_data <- separate(all_data, id, into = c("id","bin"), sep="_(?=[^_]+$)" )
all_data <- separate(all_data, id, into = c("id","exon_intron_id"), sep="_(?=[^_]+$)" )

# Remove weird extra row?
all_data <- all_data[-641,]

# Add back plotting column
all_data$plotting <- paste(all_data$id, all_data$exon_intron_id, all_data$bin, sep="_")

# Change order of features so it makes sense
all_data$bin <- as.numeric(all_data$bin)
all_data$weightedMeth.mean <- as.numeric(all_data$weightedMeth.mean)

all_data$order <- 1
all_data$order[all_data$plotting %like% "five"] <- 2
all_data$order[all_data$plotting %like% "exon_1"] <- 3
all_data$order[all_data$plotting %like% "intron_1"] <- 4
all_data$order[all_data$plotting %like% "exon_2"] <- 5
all_data$order[all_data$plotting %like% "intron_2"] <- 6
all_data$order[all_data$plotting %like% "exon_3"] <- 7
all_data$order[all_data$plotting %like% "intron_3"] <- 8
all_data$order[all_data$plotting %like% "exon_4"] <- 9
all_data$order[all_data$plotting %like% "intron_4"] <- 10
all_data$order[all_data$plotting %like% "exon_5"] <- 11
all_data$order[all_data$plotting %like% "intron_5"] <- 12
all_data$order[all_data$plotting %like% "exon_n"] <- 13
all_data$order[all_data$plotting %like% "intron_n"] <- 14
all_data$order[all_data$plotting %like% "three"] <- 15
all_data$order[all_data$plotting %like% "down"] <- 16

all_data <- all_data[order(all_data$order, all_data$bin),] 
head(all_data)

ggplot(all_data, aes(x = plotting, y= weightedMeth.mean, colour=sex))+
  geom_point()+
  geom_smooth()+
  scale_y_continuous(n.breaks = 10)+
  scale_x_discrete(labels=NULL)
  
