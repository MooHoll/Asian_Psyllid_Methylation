library(methylKit)
library(readr)

file.list <- list("F1-1_BDLM200001629-1A.bam", "F1-2_BDLM200001630-1A.bam",
                  "F1-3_BDLM200001631-1A.bam", "M1-1_BDLM200001632-1A.bam",
                  "M1-2_BDLM200001633-1A.bam", "M1-3_BDLM200001634-1A.bam")

raw_data <- processBismarkAln(file.list,
                              sample.id = list("F1","F2","F3","M1","M2","M3"),
                              treatment = c(0,0,0,1,1,1),
                              assembly="Dcitriv3", 
                              read.context="CpG",
                              save.context = "CpG",
                              save.folder=getwd(),
                              save.db = TRUE)
