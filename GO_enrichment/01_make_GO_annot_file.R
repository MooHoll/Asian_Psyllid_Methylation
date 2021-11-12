# Make a useable file from the excel of the eggNOG output

# NOTE: First removed all info from the .xlxs file except gene id and GO annotation

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/GO_enrichment")

library(tidyr)
library(dplyr)
library(readxl)

out_emapper_annotations_copy <- read_excel("Dropbox/Edinburgh/Projects/Asian_psyllid/GO_enrichment/eggNOG_GO_terms/out.emapper.annotations copy.xlsx", 
                                           col_names = FALSE)
colnames(out_emapper_annotations_copy) <- c("gene_id","GO_term")

new <- out_emapper_annotations_copy %>% 
  mutate(GO_term = strsplit(as.character(GO_term), ",")) %>% 
  unnest(GO_term)

write.table(new, file="d_citri_GO_annotations.txt", sep="\t", quote = F, row.names = F,
            col.names = T)
