## -------------------------------------------------------------------------
# Weighted methylation per window for each sex
## -------------------------------------------------------------------------

library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

# Read in sample methylation count files edited from the bismark files
file.list = list.files(("./"),pattern="*final_coverage.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

# Merge the samples from this point to save generating loads of files and merging later
females <- bind_rows(samples[1:6])
males <- bind_rows(samples[7:12])

females <- summaryBy(total_coverage + count_c ~ chr + cpg, data=females, FUN=sum)
males <- summaryBy(total_coverage + count_c ~ chr + cpg, data=males, FUN=sum)

colnames(females) <- c("chr","cpg","total_coverage","count_c")
colnames(males) <- c("chr","cpg","total_coverage","count_c")

females$chr <- as.factor(females$chr)
males$chr <- as.factor(males$chr)

females$cpg <- as.numeric(females$cpg)
males$cpg <- as.numeric(males$cpg)

both <- list(females,males)
sample_names <- list("female","male")
names(both) <- sample_names

# Read in window with start/end and total CpGs per gene
annotation_with_total_cpgs <- read_table2("windows_with_total_cpgs.txt")
colnames(annotation_with_total_cpgs) <- c("id","chr","start","end","cpg_count")
annotation_with_total_cpgs$chr <- as.factor(annotation_with_total_cpgs$chr)
annotation_with_total_cpgs$start <- as.numeric(annotation_with_total_cpgs$start)
annotation_with_total_cpgs$end <- as.numeric(annotation_with_total_cpgs$end)
## -------------------------------------------------------------------------
registerDoParallel(cores = 2)

# Calculate weighted meth for each gene for each sample
foreach(i = seq_along(both)) %dopar% {
  df <- both[[i]]
  df <- subset(df, total_coverage > 5)
  output <- sqldf("SELECT sample.chr,
                    sample.cpg,
                    sample.count_c,
                    sample.total_coverage,
                    annot.chr,
                    annot.start,
                    annot.end,
                    annot.id,
                    annot.cpg_count
                    FROM df AS sample
                    LEFT JOIN annotation_with_total_cpgs AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg >= annot.start AND sample.cpg <= annot.end)")
  output <- output[!is.na(output$id),]
  output <- output[,-c(1,2)]
  check <- summaryBy(total_coverage + count_c ~ chr + id + start + end + cpg_count, data=output, FUN=sum) 
  check$weightedMeth <- (check$cpg_count*check$count_c.sum)/(check$cpg_count*check$total_coverage.sum)
  myfile <- file.path("./", paste0(names(both[i]),"_","weighted_meth.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}


## -------------------------------------------------------------------------
