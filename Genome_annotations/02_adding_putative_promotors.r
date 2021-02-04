# ----------------------------------------------------------------
### Adding putative promotor information to the annotation file
# ----------------------------------------------------------------
setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files")
library(readr)

Dcitr_OGSv3_0_beta <- read_delim("Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
head(Dcitr_OGSv3_0_beta)
Dcitr_OGSv3_0_beta <- Dcitr_OGSv3_0_beta[,-c(2,6,8)]
Dcitr_OGSv3_0_beta <- Dcitr_OGSv3_0_beta[!duplicated(Dcitr_OGSv3_0_beta),] #0
colnames(Dcitr_OGSv3_0_beta) <- c("chr","feature","start", "end","strand","gene_id")

upsteam5UTRs <- subset(Dcitr_OGSv3_0_beta, feature == "five_prime_UTR")
upsteam5UTRs <- subset(upsteam5UTRs, strand == "+")
upsteam3UTRs <- subset(Dcitr_OGSv3_0_beta, feature == "three_prime_UTR")
upsteam3UTRs <- subset(upsteam3UTRs, strand == "-")

upsteam5UTRs$promotor_start <- upsteam5UTRs$start - 500
upsteam5UTRs <- upsteam5UTRs[upsteam5UTRs$promotor_start > 0,] #0 genes (promotor overlaps scaffold start: removed)

upsteam3UTRs$promotor_start <- upsteam3UTRs$start - 500
upsteam3UTRs <- upsteam3UTRs[upsteam3UTRs$promotor_start > 0,]

promoters <- rbind(upsteam3UTRs,upsteam5UTRs)

promoters <- promoters[,-4] # rm redundent end column
colnames(promoters)[3] <- "end"
colnames(promoters)[6] <- "start"
promoters$feature <- "promoter"

all <- rbind(Dcitr_OGSv3_0_beta,promoters)

write.table(all, file="Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters.txt ", 
            row.names = F, col.names = T, sep = '\t', quote = F)

# ----
# Add in the TE annotation as well
TEs <- read_delim("Diaci_v3.0.ref.fa.mod.EDTA.intact.txt", 
                  "\t", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE)

head(TEs)
TEs <- TEs[,-c(2,6,8)]
TEs <- TEs[!duplicated(TEs),] #2
TEs$X3 <- "TE"
colnames(TEs) <- c("chr","feature","start", "end","strand","gene_id")

all2 <- rbind(all, TEs)
write.table(all2, file="Dcitr_OGSv3.0_beta_longestIsoform_plusIntrons_plusPromoters_plusTEs.txt ", 
            row.names = F, col.names = T, sep = '\t', quote = F)
