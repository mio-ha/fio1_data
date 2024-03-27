library(edgeR)
count <- read.table("../count/dex_counts_down_flat_6.txt", sep = "\t", header = T, skip =1)
counts <- as.matrix(count[,7:10])
colnames(counts) = c("Col0.1", "Col0.2", "fio1.1", "fio1.2")
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
write.table(table, "edgeR_result_splicing_Parker_667_2rep_6.txt", col.names = T, row.names = F, sep = "\t")
