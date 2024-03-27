library(sva)
library(stringr)

#type <- c("A3SS", "A5SS", "MXE", "RI", "SE")
type <- "SE"

for (i in type){
readfile <- str_interp("rmats_all/JCEC.raw.input.${i}.txt")
outfile <- str_interp("rmats_all/${i}_combat_adj.txt")

A3 <- read.table(text = gsub(",", "\t", readLines(readfile)), skip=1, fill=T)
A3c <- A3[, c(2:57)]

batch <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3)
cov1 <- as.numeric(factor(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
cov2 <- as.numeric(factor(c(rep(0, 13), rep(1, 13), rep(0, 15), rep(1, 15))))
cov_m <- model.matrix(as.formula(~0+cov1+cov2), data=A3c)

A3adj <- ComBat(dat=A3c, batch=batch, mod=cov_m, par.prior=TRUE, mean.only=TRUE)
A3merge <- cbind(A3[, 1], A3adj, A3[, c(58,59)])

A3out <- as.data.frame(lapply(A3merge, function(x) as.integer(x)))
A3out[A3out < 0] <- 0

write.table(A3out, file=outfile, sep="\t", row.names = F, col.names = F)
}
