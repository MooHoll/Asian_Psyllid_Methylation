## -------------------------------------------------------------------------
## Weighted Methylation per window: make file with cpg count per window
## -------------------------------------------------------------------------

# NOTE: run this on HPC, huge memory and processing needed

# Load packages etc.
library(sqldf)
library(readr)
library(doBy)
library(dplyr)
library(foreach)
library(doParallel)

## -------------------------------------------------------------------------
# Making the file which has the total CpGs per gene information

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("chr", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)
cpgs$chr <- as.factor(cpgs$chr)

# --------------------------------------------------------------------
# List files starting with xa (from the split command)
file.list = list.files(("./"),pattern="^xa*")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

samples <- lapply(file.list, read_file1)

for(i in seq_along(samples)){
  colnames(samples[[i]]) <- c("chr","start","end","id")
  samples[[i]]$chr <- as.factor(samples[[i]]$chr)
  samples[[i]]$start <- as.numeric(samples[[i]]$start)
  samples[[i]]$end <- as.numeric(samples[[i]]$end)
}

registerDoParallel(cores = 20)

list_of_outputs <- list()

list_of_outputs <- foreach(i = seq_along(samples)) %dopar% {
  windows <- samples[[i]]
  output <- sqldf("SELECT cpgs.chr,
                cpgs.cpg_position,
                windows.chr,
                windows.start,
                windows.end,
                windows.id
                FROM windows AS windows
                LEFT JOIN cpgs AS cpgs 
                ON cpgs.chr = windows.chr
                AND (cpgs.cpg_position >= windows.start AND cpgs.cpg_position <= windows.end)")
  return(output)
}

final <- bind_rows(list_of_outputs)
final <- final[!is.na(final$id),]
final$cpg_counter <- 1
final <- final[,-c(1:2)]
colnames(final) <- c("chr", "start","end","id","cpg_counter")

final <- summaryBy(cpg_counter ~ id+chr+start+end, data = final, FUN=sum)

write.table(final, file="windows_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)
