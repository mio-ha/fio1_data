library(edgeR)

for (i in 1:6) {
for (j in 1:6) {

if (i >= j) {
	next
}

open.file <- paste0("../count/2rep_matrix/dex_counts_Parker_2rep_", i, "_", j, ".txt")
print(open.file)

count <- read.table(open.file, sep = "\t", header = T, skip =1)
counts <- as.matrix(count[,7:10])
colnames(counts) = c(paste0("Col-0_", i), paste0("Col-0_", j), paste0("fio1_", i), paste0("fio1_", j))
genes <- count[,1:6]

y <- DGEList(counts=counts, genes=genes)
group = factor(c("col0", "col0", "fio1", "fio1"))
y$samples$group <- group

keep <- filterByExpr(y, group=group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~ 0 + group)

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)

contr <- makeContrasts(groupfio1 - groupcol0, levels=design)
sp <- diffSpliceDGE(fit, contrast=contr, geneid="Geneid", exonid="Start")
table <- as.data.frame(topSpliceDGE(sp, test="exon", n=nrow(sp$genes)))
write.table(table, paste0("2rep/edgeR_result_splicing_Parker_2rep_", i, "_", j, ".txt"), col.names = T, row.names = F, sep = "\t")

}
}
