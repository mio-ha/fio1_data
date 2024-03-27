library(sva)
library(stringr)

type <- c("A3SS", "A5SS", "MXE", "RI", "SE")

for (i in type){
readfile <- str_interp("rmats_all/JCEC.raw.input.${i}.txt")
outfile <- str_interp("rmats_all/${i}_combat_adj.txt")

A3 <- read.table(text = gsub(",", "\t", readLines(readfile)), skip=1, fill=T)
A3c <- A3[, c(2:57)]

batch <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3)
cov1 <- as.numeric(factor(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
cov2 <- as.numeric(factor(c(rep(0, 13), rep(1, 13), rep(0, 15), rep(1, 15))))
covar_mat <- cbind(cov1, cov2)

A3adj <- ComBat_seq(A3c, batch=batch, group=NULL, covar_mod=covar_mat)
A3out <- cbind(A3[, 1], A3adj, A3[, c(58,59)])
write.table(A3out, file=outfile, sep="\t", row.names = F, col.names = F)
}
