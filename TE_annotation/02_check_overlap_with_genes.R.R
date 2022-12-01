# Check that TE annotations don't significantly overlap with gene annotations

setwd("~/Dropbox/Edinburgh/Projects/Asian_psyllid/Genome_Files")

library(readr)
library(sqldf)
library(reshape2)
library(ggplot2)

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

# Should probably just confirm as well that the methylation levels between
# gene TEs and non-gene TEs are similar
Dcitri_weighted_meth_annotation_by_sex <- read_delim("~/Dropbox/Edinburgh/Projects/Asian_psyllid/DNA_Methylation/Differential_methylation/Dcitri_weighted_meth_annotation_by_sex.txt", 
                                                     delim = "\t", escape_double = FALSE, 
                                                     trim_ws = TRUE)
TE_meth <- Dcitri_weighted_meth_annotation_by_sex[Dcitri_weighted_meth_annotation_by_sex$feature=="TE",]
head(TE_meth)
TE_meth <- TE_meth[,-c(2,6)]


both <- rbind(output_start, output_end)
head(both)
colnames(both)[10]<-"actual_gene_id"
both$in_gene <-"no"
both$in_gene[!is.na(both$actual_gene_id)]<-"yes"
table(both$in_gene)
both <- both[,-c(2,6,7,8,9,10,11)]

both <- both[!duplicated(both),] 
both$unique <- paste(both$start, both$end)
both <- both[!duplicated(both$unique),] 
table(both$in_gene)
#no  yes 
#6068 6093 

look <- merge(both, TE_meth, by=c("chr","start",'end'))
look <- look[,c(4,5,8,9)]
look_long <- melt(look, id.vars = c("gene_id.x","in_gene"))
colnames(look_long) <- c("TE","in_gene","Sex","meth")
look_long$log_meth <- log10(look_long$meth)


ggplot(look_long, aes(x=in_gene, y=log_meth))+
  geom_violin(aes(fill=Sex))+
  xlab("Location")+
  ylab("log(Weighted Methylation)")+
  scale_fill_manual(breaks = c("female","male"),labels=c("Female","Male"),
                    values=c("#DDCC77","#44AA99"))+
  scale_x_discrete(limits=c("no","yes"),
                   labels=c("No Gene Overlap", "Within/Overlapping a Gene"))+
  theme_bw()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=20),
        axis.title.x = element_text(size=20),
        plot.title=element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank())

# Stats
library(multcomp)
head(look_long)
model1<-lm(meth ~ Sex * in_gene , data=look_long)
model2<-lm(meth ~ Sex + in_gene , data=look_long)
anova(model1,model2) # No interaction
summary.lm(model2) # In gene does have an effect on the level of methylation

mean(look_long$meth[look_long$in_gene=="no"], na.rm = T)
mean(look_long$meth[look_long$in_gene=="yes"], na.rm = T)
median(look_long$meth[look_long$in_gene=="no"], na.rm = T)
median(look_long$meth[look_long$in_gene=="yes"], na.rm = T)

