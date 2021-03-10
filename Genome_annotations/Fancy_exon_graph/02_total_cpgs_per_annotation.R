## -------------------------------------------------------------------------
## Weighted Methylation per window: make file with cpg count per window
## -------------------------------------------------------------------------

# NOTE: run this on HPC, huge memory and processing needed

# Load packages etc.
library(sqldf)
library(readr)
library(doBy)

## -------------------------------------------------------------------------
# Making the file which has the total CpGs per gene information

cpgs <- read.delim("total_cpgs_in_genome.txt", header=F)
colnames(cpgs) <- c("chr", "cpg_position")
cpgs$cpg_position <- as.numeric(cpgs$cpg_position)
head(cpgs)

# --------------------------------------------------------------------

windows <- read.delim("windows.bed", header=F)
colnames(windows) <- c("chr","start","end","id")
windows$start <- as.numeric(windows$start)
windows$end <- as.numeric(windows$end)

output <- sqldf("SELECT cpgs.chr,
                cpgs.cpg_position,
                windows.chr,
                windows.start,
                windows.end,
                windows.id
                FROM cpgs AS cpgs
                LEFT JOIN windows AS windows 
                ON cpgs.chr = windows.chr
                AND (cpgs.cpg_position >= windows.start AND cpgs.cpg_position <= windows.end)")
output <- output[!is.na(output$id),]

output$cpg_counter <- 1

final <- summaryBy(cpg_counter ~ id+chr+start+end, data = output, FUN=sum)

write.table(final, file="windows_with_total_cpgs.txt", col.names=T, row.names=F, quote=F)
