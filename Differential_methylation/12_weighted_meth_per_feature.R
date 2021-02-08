## -------------------------------------------------------------------------
# Weighted methylation per annotation for each sex
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

# Read in sample methylation count files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("F1_R1", "F1_R2", "F2_R1","F2_R2", "F3_R1", "F3_R2",
                    "M1_R1", "M1_R2", "M2_R1","M2_R2", "M3_R1", "M3_R2")
names(samples) <- sample_names

# Read in gene with start/end and total CpGs per gene
annotation_with_total_cpgs <- read_table2("Dcitr_OGSv3.0_beta_with_total_cpgs.txt")

## -------------------------------------------------------------------------
registerDoParallel(cores = 10)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.cpg_count,
                    annot.feature
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + feature + gene_id + start + end + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(samples[i]),"_","weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}


## -------------------------------------------------------------------------
