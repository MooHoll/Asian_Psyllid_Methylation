## -------------------------------------------------------------------------
# Making scatter with diff meth genes highlighted
## -------------------------------------------------------------------------
setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation")
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)

all_with_meth <- read_delim("Dcitri_weighted_meth_annotation_by_sex.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

diff_meth_exon_geneIDs <- read_csv("diff_meth_exon_geneIDs.txt", 
                                   col_names = FALSE)
colnames(diff_meth_exon_geneIDs) <- "gene_id"

## -------------------------------------------------------------------------

exon_data <- all_with_meth[all_with_meth$feature=="exon" ,]

exon_data$diff <- "no"
exon_data$diff[exon_data$gene_id %in% diff_meth_exon_geneIDs$gene_id] <- "yes"
write.table(exon_data, file="all_exon_information.txt", sep="\t",
            quote = F, col.names = T, row.names = F)


ggplot(exon_data %>%
         arrange(diff), aes(x=female, y=male))+
  geom_point(aes(colour=diff), size=2.25)+
  ylim(0,1.0)+
  theme_bw()+
  xlab("Female Weighted Methylation")+
  ylab("Male Weighted Methylation")+
  ggtitle("Exons")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title=element_text(size = 22),
        legend.position = "none")+
  scale_colour_manual(breaks = c("no","yes"),
                      values=c("black","red"))

