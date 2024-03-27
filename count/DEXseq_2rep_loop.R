source("/home/miyokawa/ドキュメント/ngs/software/Subread_to_DEXSeq/load_SubreadOutput.R")
suppressPackageStartupMessages({
  require(dplyr)
})

for (i in 1:6) {
for (j in 1:6) {

if (i >= j) {
	next
}

# Load featureCounts result file
rn <- c(paste0("Parker_Col-0_", i, ".bam"), paste0("Parker_Col-0_", j, ".bam"), paste0("Parker_fio1_", i, ".bam"), paste0("Parker_fio1_", j, ".bam"))
samp <- data.frame(row.names = rn, 
                   condition = c("Col_0", "Col_0", "fio1", "fio1"))
open.file <- paste0("2rep_matrix/dex_counts_Parker_2rep_", i, "_", j, ".txt")
print(open.file)
dxd.fc <- DEXSeqDataSetFromFeatureCounts(open.file,
                                         flattenedfile = "Araport11_GTF_genes_transposons.Apr2023_flat.gtf",
                                         sampleData = samp)
# Adjust for coverage biases among samples
dxd <- estimateSizeFactors(dxd.fc)

# Estimate variance or dispersion parameters individually exon by exon
dxd <- estimateDispersions(dxd)

#png(filename = "test.png", width = 600, height = 600)
#plotDispEsts(dxd)
#dev.off()

# Test for differential exon usage in each gene
dxd <- testForDEU(dxd)

# Estimate relative exon usage fold changes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")

# Output DEXSeq result
dxr1 <- DEXSeqResults(dxd)
table <- as.data.frame(dxr1)
table$transcripts <- as.character(table$transcripts)
write.table(table, quote = F, sep = "\t", file = paste0("../DEXSeq/DEXSeq_Parker_2rep_", i, "_", j, "_result.txt"), row.names = F)

}
}
