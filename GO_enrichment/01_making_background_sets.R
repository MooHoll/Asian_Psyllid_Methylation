# ------------------------------------------------------
# Making base genome-wide GO file and background sets
# ------------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/GO_enrichment")
library(readr)
library(tidyr)

# ------------------------------------------------------
# Make a compatible GO annotation file (only 5861 annotated GO terms)
annotations <- read_delim("Dcitr_OGSv3.0_beta_pep.desc.fa.goanna.inverteb.exponly.percentident70_qcov70.out_goanna_gaf.tsv", 
                          "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
annotations <- annotations[,c(2,5)] 
colnames(annotations) <- c("geneID","goID")
annotations$geneID <- gsub("\\.1.*", ".1", annotations$geneID)
new_annotations <- separate_rows(annotations, goID, sep =';')
write.table(new_annotations, file="D_citri_GO_terms.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

# ------------------------------------------------------
# Make gene lists

# Differentially expressed genes
all_gene_expression_data <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Gene_Expression/Differential_expression/all_gene_expression_data.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)

all_RNA_genes <- as.data.frame(all_gene_expression_data$gene_id) # 12420
colnames(all_RNA_genes) <- "geneID"

diff_exp_genes <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$diff_exp=="yes"]) #1259
colnames(diff_exp_genes) <- "geneID" 

male_biased <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category=="male_biased"]) # 675
colnames(male_biased)<- "geneID"

male_limited <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category=="male_limited"]) # 484
colnames(male_limited)<- "geneID"

male_biased_extreme <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category=="male_biased_extreme"]) # 5
colnames(male_biased_extreme)<- "geneID"

all_male_biased_categories <- rbind(male_biased, male_biased_extreme, male_limited) # 1164

female_biased <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category=="female_biased"]) # 83
colnames(female_biased)<- "geneID"

female_limited <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category=="female_limited"]) # 12
colnames(female_limited)<- "geneID"

all_female_biased_categories <- rbind(female_biased, female_limited) #95

# ------------------------------------------------------
# Make background sets

all_RNA_genes_background <- merge(annotations, all_RNA_genes) 
length(unique(all_RNA_genes_background$geneID)) # 338/12420
write.table(as.data.frame(all_RNA_genes_background), file="./background_go_sets/all_RNA_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

all_diff_exp_genes_background <- merge(annotations, diff_exp_genes)
length(unique(all_diff_exp_genes_background$geneID)) # 10/1259
write.table(as.data.frame(all_diff_exp_genes_background), file="./background_go_sets/all_diff_exp_genes_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

all_male_biased_categories_background <- merge(annotations, all_male_biased_categories)
length(unique(all_male_biased_categories_background$geneID)) # 10/1164
write.table(as.data.frame(all_male_biased_categories_background), file="./background_go_sets/all_male_biased_categories_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)

all_female_biased_categories_background <- merge(annotations, all_female_biased_categories) 
length(unique(all_female_biased_categories_background$geneID))# 0/95
write.table(as.data.frame(all_female_biased_categories_background), file="./background_go_sets/all_female_biased_categories_background.txt", sep="\t", quote = F,
            col.names = T, row.names = F)
