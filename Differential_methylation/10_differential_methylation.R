## -------------------------------------------------------------------------
## Differential Methylation Between Female and Male D.citri
## -------------------------------------------------------------------------

# Load packages etc.
library(methylKit)
library(readr)

## -------------------------------------------------------------------------
# Put in the correct comparison below, see script: 08_comparisons_diff_meth.R

file.list <- list("F1_CpG.txt", "F2_CpG.txt",
                  "F3_CpG.txt", "M1_CpG.txt",
                  "M2_CpG.txt", "M3_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("F1","F2","F3","M1","M2","M3"),
                     treatment = c(0,0,0,1,1,1),
                     assembly="Dcitriv3", 
                     context="CpG")

## -------------------------------------------------------------------------

# Filter data for outliers and coverage and normalise counts
filtered_data <- filterByCoverage(raw_data,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

normalized <- normalizeCoverage(filtered_data)

## -------------------------------------------------------------------------

# Only use CpGs present in all samples (maybe too stringent, can adjust later)
meth_all_data <- unite(normalized, destrand=TRUE) 
nrow(meth_all_data) # 3206305

## -------------------------------------------------------------------------

# Filter sites using a binomial test so only keep CpGs which are methylated in at least one sample
df_meth_all <- getData(meth_all_data)

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]
c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]
e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]

# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.006) {binom.test(a, b, 0.006, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 107710

subset_methBase <- methylKit::select(meth_all_data, meth_positions)

## -------------------------------------------------------------------------

# Save the dataframe for later use, including making nice plots
methBase_ob <- getData(subset_methBase)
write.table(methBase_ob, file="F_vs_M_objectmethbase.txt", quote=F, row.names = F, sep = '\t')

## -------------------------------------------------------------------------

# Generate some quick plots which can be fixed up later
pdf("F_vs_M_correlation.pdf")
getCorrelation(subset_methBase,plot=TRUE)
dev.off()

pdf("F_vs_M_cluster.pdf")
clusterSamples(subset_methBase, dist="correlation", method="ward", plot=TRUE)
dev.off()

pdf("F_vs_M_PCA.pdf")
PCASamples(subset_methBase)
dev.off()


## -------------------------------------------------------------------------

# Differential methylation as group 1 minus group 0
diff_meth <- calculateDiffMeth(subset_methBase, mc.cores = 1)
write.csv(diff_meth, file="F_vs_M__all_tested_meth_sites_MSCfilter.csv")

diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.05)
write.csv(diff_meth_10, file="F_vs_M__DMRs_min10percentDiff_qval0.05.csv")
nrow(diff_meth_10) # 3112/107710

diff_meth_10 <- getMethylDiff(diff_meth, difference=10, qvalue=0.01)
write.csv(diff_meth_10, file="F_vs_M__DMRs_min10percentDiff_qval0.01.csv")
nrow(diff_meth_10) # 763/107710

diff_meth_15 <- getMethylDiff(diff_meth, difference=25, qvalue=0.01)
write.csv(diff_meth_15, file="F_vs_M__DMRs_min25percentDiff_qval0.01.csv")
nrow(diff_meth_15) # 261/107710

# Use this threshold
diff_meth_15 <- getMethylDiff(diff_meth, difference=15, qvalue=0.01)
write.csv(diff_meth_15, file="F_vs_M__DMRs_min15percentDiff_qval0.01.csv")
nrow(diff_meth_15) # 600/107710

