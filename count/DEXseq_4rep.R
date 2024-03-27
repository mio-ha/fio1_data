source("/home/miyokawa/ドキュメント/ngs/software/Subread_to_DEXSeq/load_SubreadOutput.R")
suppressPackageStartupMessages({
  require(dplyr)
})

# Load featureCounts result file
rn <- c("Parker_Col-0_1.bam","Parker_Col-0_2.bam", "Parker_Col-0_3.bam", "Parker_Col-0_6.bam", "Parker_fio1_1.bam", "Parker_fio1_2.bam", "Parker_fio1_4.bam", "Parker_fio1_5.bam")
samp <- data.frame(row.names = rn, 
                   condition = c("Col_0", "Col_0", "Col_0", "Col_0", "fio1", "fio1", "fio1", "fio1"))
dxd.fc <- DEXSeqDataSetFromFeatureCounts("dex_counts_Parker20c_4rep.txt",
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
write.table(table, quote = F, sep = "\t", file = "../DEXSeq/DEXSeq_Parker20c_4rep_result.txt", row.names = F)
