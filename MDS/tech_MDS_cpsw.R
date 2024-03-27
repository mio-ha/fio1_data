library(edgeR)

count <- read.table("count/dex_counts_cpsw.csv", sep = "\t", header = T, row.names = 1)
count <- as.matrix(count)

colnames(count) = c("Cai_Col-0_1", "Cai_Col-0_2", "Cai_Col-0_3", "Cai_fio1_1", "Cai_fio1_2", "Cai_fio1_3", "Parker_Col-0_20c_1", "Parker_Col-0_20c_2", "Parker_Col-0_20c_3", "Parker_Col-0_20c_4", "Parker_Col-0_20c_5", "Parker_Col-0_20c_6", "Parker_fio1_20c_1", "Parker_fio1_20c_2", "Parker_fio1_20c_3", "Parker_fio1_20c_4", "Parker_fio1_20c_5", "Parker_fio1_20c_6", "Sun_Col-0_1", "Sun_Col-0_2", "Sun_fio1-1_1", "Sun_fio1-1_2", "Sun_fio1-5_1", "Sun_fio1-5_2", "Wang_Col-0_1", "Wang_Col-0_2", "Wang_fio1_1", "Wang_fio1_2")
group = factor(c("Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1"))
treat = factor(c("Cai", "Cai", "Cai", "Cai", "Cai", "Cai", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Wang", "Wang", "Wang", "Wang"))

design <- model.matrix(~ group + treat)
d <- DGEList(counts = count, group = group)
keep <- filterByExpr(d, design)
d <- d[keep,,keep.lib.sizes = F]
d <- calcNormFactors(d)
logcpm <- cpm(d, log = T)

designb <- model.matrix(~ group)
logcpm_nob <- removeBatchEffect(logcpm, batch=treat, design=designb)

color <- c("orangered", "mediumaquamarine", "skyblue", "royalblue")[treat]
marker <- c(1, 19)[group]
#print(head(logcpm))
pdf("MDS_cpsw_nobatch2.pdf")
plotMDS(logcpm_nob, col=color, pch=marker, cex=1.5)
#plotMDS(logcpm_nob, col=color, pch=marker)
legend("topright",fill=c("orangered", "mediumaquamarine", "skyblue", "royalblue"), legend=levels(treat))
legend("topleft",pch = c(1, 19),legend=levels(group))
dev.off()
