# ------------------------------------------------------
# Making that fancy exon and intron methylation graph
# ------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files")

library(readr)
library(dplyr)
library(tidyr)

# ------------------------------------------------------
# check the filke is sorted so the 1st exon is 1st etc. 
# also check for strandedness as this could mess up the order
annotation <- read_delim("Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters_plusTEs_plusIntergenic.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)

# ------------------------------------------------------
# Take only the exons
exons_only <- subset(annotation, feature=="exon")

# Number each exon per gene
numbered_exons <- exons_only %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

# Label any exons >=6 as just n
numbered_exons$number[numbered_exons$number>=6] <- "n"

# Do the same for introns
introns_only <- subset(annotation, feature=="intron")

numbered_introns <- introns_only %>%
  group_by(gene_id) %>%
  mutate(number = 1:n())

numbered_introns$number[numbered_introns$number>=6] <- "n"

# ------------------------------------------------------
# Replace the exons and introns in the main file

# Remove the original exons and introns
annotation1 <- annotation[!(annotation$feature == "exon" | annotation$feature=="intron"),]

# Put in the new ones
annotation1$number <- NA
final_annotation <- rbind(annotation1, numbered_exons, numbered_introns)

#Write out the file
write.table(final_annotation, file ="Dcitr_OGSv3.0_beta_numbered_exons.txt", sep="\t",
            quote = F, col.names = T, row.names = F)

# ------------------------------------------------------
# Get all the data together for the plot
genes <- final_annotation[final_annotation$feature=="gene",]

# again watch out incase your data is strand specific here, if so you need to do something different
genes$upstream <- genes$start - 5000

# Make a new annotation for an upstream region
upstream <- genes
upstream$feature <- "upstream_region"
upstream <- upstream[,-4]
colnames(upstream)[3]<- "end"
colnames(upstream)[6]<- "start"

# Make a new annotation for a downstream region
genes <- final_annotation[final_annotation$feature=="gene",]
genes$downstream <- genes$end + 5000
downstream <- genes
downstream$feature <- "downstream_region"
downstream <- downstream[,-3]
colnames(downstream)[3]<- "start"
colnames(downstream)[6]<- "end"

# Get all the data together that we want to plot
head(upstream)
head(downstream)
head(numbered_exons)
head(numbered_introns)

UTRs <- final_annotation[(final_annotation$feature == "five_prime_UTR" | 
                            final_annotation$feature == "three_prime_UTR"),]
head(UTRs)

plotting_data <- rbind(upstream, downstream, numbered_exons, numbered_introns, UTRs)
plotting_data$number <- replace_na(plotting_data$number, 1)

# Remove any rows where the start is <0 i.e. overlaps the start of a chromosome
plotting_data <- plotting_data[plotting_data$start > 0,]

# ------------------------------------------------------
# Write out file in bed format compatible with bedtools makewindows
plotting_data <- plotting_data[,c(1,6,3,2,5)]
plotting_data$id <- paste0(plotting_data$feature,"_", plotting_data$number)
plotting_data <- plotting_data[,c(1:3,6)]

write.table(plotting_data, file ="Dcitr_OGSv3.0_beta_for_exon_graph.txt", sep="\t",
            quote = F, col.names = F, row.names = F)

# ------------------------------------------------------
# In bash run this to make the windows:
# bedtools makewindows -b Dcitr_OGSv3.0_beta_for_exon_graph.txt -n 20 -i srcwinnum > windows.bed

