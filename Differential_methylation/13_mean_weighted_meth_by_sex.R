## -------------------------------------------------------------------------
# Take average weighted methylation level of feature across bio replicates
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(reshape2)

# Make one file covering all samples
file.list = list.files(("./"),pattern="*weighted_meth.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)


# Make one dataframe for each population
females <- samples[1:6]
males <- samples[7:12]

for(i in seq_along(females)){
  females[[i]]$sex <- "female"
}
females_all <- as.data.frame(bind_rows(females))
females_merged <- summaryBy(weightedMeth ~ chr + feature + gene_id + start + end +
                              cpg_count + sex, data = females_all, FUN=mean)

for(i in seq_along(males)){
  males[[i]]$sex <- "male"
}
males_all <- as.data.frame(bind_rows(males))
males_merged <- summaryBy(weightedMeth ~ chr + feature + gene_id + start + end +
                            cpg_count + sex, data = males_all, FUN=mean)

all_data <- rbind(females_merged, males_merged)

all_data2 <- dcast(all_data, chr + feature + gene_id + start + end +
                     cpg_count ~ sex, value.var = "weightedMeth.mean")

write.table(all_data2, file="Dcitri_weighted_meth_annotation_by_sex.txt", col.names = T,
            row.names = F, quote = F, sep = "\t")



