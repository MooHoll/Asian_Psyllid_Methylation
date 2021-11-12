#---------------------------------------------------
# Background gene sets for differential meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/GO_enrichment")
library(readr)
library(stringr)

new_d_citri_GO_annotations <- read_delim("new_d_citri_GO_annotations.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
new_d_citri_GO_annotations$gene_id <- str_sub(new_d_citri_GO_annotations$gene_id, end=-3)

#---------------------------------------------------

methylated_genes <- read_delim("background_go_sets/methylated_genes.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
colnames(methylated_genes) <- "gene_id"
methylated_genes <- methylated_genes[!duplicated(methylated_genes),]

meth_genes_with_go <- merge(methylated_genes, new_d_citri_GO_annotations, by = "gene_id")
meth_genes_with_go <- meth_genes_with_go[!meth_genes_with_go$GO_term=="-",]
length(unique(meth_genes_with_go$gene_id)) #3025/4429

write.table(meth_genes_with_go, file="./background_go_sets/methylated_genes_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------

all_gene_expression_data <- read_delim("background_go_sets/all_gene_expression_data.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
head(all_gene_expression_data)

all_gene_expression_data <- all_gene_expression_data[,c(1,16,17)]
diff_exp <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$diff_exp =="yes"])
colnames(diff_exp) <- "gene_id"
diff_exp_genes_with_go <- merge(diff_exp, new_d_citri_GO_annotations, by = "gene_id")
diff_exp_genes_with_go <- diff_exp_genes_with_go[!diff_exp_genes_with_go$GO_term=="-",]
length(unique(diff_exp_genes_with_go$gene_id)) #462/1259
write.table(diff_exp_genes_with_go, file="./background_go_sets/diff_exp_genes_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

diff_exp_list <- as.data.frame(unique(diff_exp$gene_id))
colnames(diff_exp_list) <- "gene_id"
write.table(diff_exp_list, file="./gene_lists/diff_exp_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)

all_rna_genes <- as.data.frame(all_gene_expression_data$gene_id[!duplicated(all_gene_expression_data$gene_id)])
colnames(all_rna_genes) <- "gene_id"
all_rna_genes_with_go <- merge(all_rna_genes, new_d_citri_GO_annotations, by = "gene_id")
all_rna_genes_with_go <- all_rna_genes_with_go[!all_rna_genes_with_go$GO_term=="-",]
length(unique(all_rna_genes_with_go$gene_id)) # 7850/12420
write.table(all_rna_genes_with_go, file="./background_go_sets/all_rna_genes_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
## Pull out the diff exp gene lists

female_biased <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category =="female_biased" |
                                                                  all_gene_expression_data$category =="female_limited"])
colnames(female_biased) <- "gene_id"
write.table(female_biased, file="./gene_lists/diff_exp_female_biased.txt", sep="\t",quote = F, col.names = T, row.names = F)


male_biased <- as.data.frame(all_gene_expression_data$gene_id[all_gene_expression_data$category =="male_biased" |
                                                                all_gene_expression_data$category =="male_limited" |
                                                                all_gene_expression_data$category == "male_biased_extreme"])
colnames(male_biased) <- "gene_id"
write.table(male_biased, file="./gene_lists/diff_exp_male_biased.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
## Pull out the diff meth gene lists

hypermethylated_genes_with_category <- read_delim("gene_lists/hypermethylated_genes_with_category.txt", 
                                                      delim = "\t", escape_double = FALSE, 
                                                      trim_ws = TRUE)
head(hypermethylated_genes_with_category) #14

diff_meth <- as.data.frame(unique(hypermethylated_genes_with_category$gene_id))
colnames(diff_meth) <- "gene_id"
diff_meth_genes_with_go <- merge(diff_meth, new_d_citri_GO_annotations, by = "gene_id")
diff_meth_genes_with_go <- diff_meth_genes_with_go[!diff_meth_genes_with_go$GO_term=="-",]
length(unique(diff_meth_genes_with_go$gene_id)) # 9/12
write.table(diff_meth_genes_with_go, file="./background_go_sets/diff_meth_genes_with_GO.txt", sep="\t",quote = F, col.names = T, row.names = F)

diff_meth_list <- as.data.frame(unique(diff_meth$gene_id))
colnames(diff_meth_list) <- "gene_id"
write.table(diff_meth_list, file="./gene_lists/diff_meth_genes.txt", sep="\t",quote = F, col.names = T, row.names = F)


hyper_male <- as.data.frame(hypermethylated_genes_with_category$gene_id[hypermethylated_genes_with_category$hypermethylated=="male"])
colnames(hyper_male) <- "gene_id"
write.table(hyper_male, file="./gene_lists/hyper_meth_male.txt", sep="\t",quote = F, col.names = T, row.names = F)

hyper_female <- as.data.frame(hypermethylated_genes_with_category$gene_id[hypermethylated_genes_with_category$hypermethylated=="female"])
colnames(hyper_female) <- "gene_id"
hyper_female <- as.data.frame(hyper_female[!duplicated(hyper_female$gene_id),])
colnames(hyper_female) <- "gene_id"
write.table(hyper_female, file="./gene_lists/hyper_meth_female.txt", sep="\t",quote = F, col.names = T, row.names = F)

