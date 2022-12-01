# Check that TE annotations don't significantly overlap with gene annotations

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files")

library(readr)
library(sqldf)

Dcitr_OGSv3_0_beta_numbered_exons <- read_delim("Dcitr_OGSv3.0_beta_numbered_exons.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)

TEs <- subset(Dcitr_OGSv3_0_beta_numbered_exons, feature=="TE")
genes <- subset(Dcitr_OGSv3_0_beta_numbered_exons, feature=="gene")

# For the start coordinate
output_start <- sqldf("SELECT sample.chr,
                    sample.feature,
                    sample.start,
                    sample.end,
                    sample.gene_id,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM TEs AS sample
                    LEFT JOIN genes AS annot
                    ON sample.chr = annot.chr
                    AND (sample.start >= annot.start AND sample.start <= annot.end)")
output_start$where <- "start"

# For the end coordinate
output_end <- sqldf("SELECT sample.chr,
                    sample.feature,
                    sample.start,
                    sample.end,
                    sample.gene_id,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id
                    FROM TEs AS sample
                    LEFT JOIN genes AS annot
                    ON sample.chr = annot.chr
                    AND (sample.end >= annot.start AND sample.end <= annot.end)")
output_end$where <- "end"

both <- rbind(output_start, output_end)

# Figure out real overlap, i.e. remove doubles caused by start and end
both <- both[,c(3,4,5,6,7,10)]
both <- subset(both, feature=="gene")
both <- both[!duplicated(both),] 

both$unique <- paste(both$start, both$end)
both <- both[!duplicated(both$unique),] 

# Hypergeometric test for overlap
# 6268 TEs and genes overlap
# 5897 TEs without gene overlap
# 12165 individual annotated TEs
# 19049 individual annotated genes

# Can say the percentage of TEs which have some level of gene overlap
100*(6268/12165) #51.52

# Next look at the total sequence percentage overlap
both <- rbind(output_start, output_end)
both <- both[,-c(2,11)]
both <- subset(both, feature=="gene")
both <- both[!duplicated(both),] 

both$unique <- paste(both$start, both$end)
both <- both[!duplicated(both$unique),] 

just_coords <- both[c(2,3,7,8)]
colnames(just_coords) <- c("TE_start","TE_end","gene_start","gene_end")

just_coords$TE_length <- just_coords$TE_end - just_coords$TE_start
range(just_coords$TE_length)
sum(just_coords$TE_length) # 8266551

just_coords$gene_length <- just_coords$gene_end - just_coords$gene_start
range(just_coords$gene_length)
sum(just_coords$gene_length) # 398077222 (this probably isn't that relevant)

# Need to figure out the length of each TE which overlaps a gene

# is the TE fully within a gene or just overlapping?
just_coords$position <- "overlapping"
just_coords$position[just_coords$TE_start > just_coords$gene_start &
                       just_coords$TE_end < just_coords$gene_end] <- "within"
table(just_coords$position) 
#overlapping      within 
# 411        5857 

overlapping <- subset(just_coords, position =="overlapping")
within <- subset(just_coords, position =="within")
sum(within$TE_length) # 6397792

overlapping$length1 <- overlapping$gene_end - overlapping$TE_start 
overlapping$length1[overlapping$length1 >= overlapping$TE_length] <- NA
overlapping$length2 <- overlapping$TE_end - overlapping$gene_start 
overlapping$length2[overlapping$length2 >= overlapping$TE_length] <- NA

sum(overlapping$length1,na.rm = T) # 435445
sum(overlapping$length2,na.rm = T) # 391608

435445 + 391608 + 6397792 # 7224845

100*(7224845/8266551) # 87% of TE of TE length of TEs which do overlap with genes

TEs$length <- TEs$end - TEs$start
sum(TEs$length) # 15658460
100*(7224845/15658460) # 46% of all TE length overlaps with genic regions

# Need to check out specifically what's going on with copia
copia <- subset(TEs, gene_id =="LTR/Copia") #186

both <- rbind(output_start, output_end)
head(both)
both <- both[,-c(2,11)]
both <- subset(both, feature=="gene")
both <- both[!duplicated(both),] 
both$unique <- paste(both$start, both$end)
both <- both[!duplicated(both$unique),] 

both_copia <- subset(both, gene_id =="LTR/Copia") # 79

just_coords <- both_copia[c(2,3,7,8)]
colnames(just_coords) <- c("TE_start","TE_end","gene_start","gene_end")

just_coords$TE_length <- just_coords$TE_end - just_coords$TE_start
range(just_coords$TE_length)
sum(just_coords$TE_length) # 136006

just_coords$position <- "overlapping"
just_coords$position[just_coords$TE_start > just_coords$gene_start &
                       just_coords$TE_end < just_coords$gene_end] <- "within"
table(just_coords$position) 
#overlapping      within 
# 11        68 

overlapping <- subset(just_coords, position =="overlapping")
within <- subset(just_coords, position =="within")
sum(within$TE_length) # 101215

overlapping$length1 <- overlapping$gene_end - overlapping$TE_start 
overlapping$length1[overlapping$length1 >= overlapping$TE_length] <- NA
overlapping$length2 <- overlapping$TE_end - overlapping$gene_start 
overlapping$length2[overlapping$length2 >= overlapping$TE_length] <- NA

sum(overlapping$length1,na.rm = T) # 22801
sum(overlapping$length2,na.rm = T) # 8228

101215 + 22801 + 8228 # 132244

copia$length <- copia$end - copia$start
sum(copia$length) # 314028
100*(132244/314028) #42.11

# Goodness of fit
observed = c(42.11, 46.14)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected)

