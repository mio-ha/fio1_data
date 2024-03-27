library(limma)
readfile <- "all_tpm_20c.txt"
outfile <- "all_tpm_20c_adj.txt"

df <- read.table(readfile, sep="\t", row.names=1, header=T)
df <- as.matrix(df)

batch <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3)
cov <- as.numeric(factor(c(rep(0, 13), rep(0, 15))))
dfadj <- removeBatchEffect(df, batch=batch, covariates=cov)
write.table(dfadj, file=outfile, sep="\t", row.names = T, col.names = T)
