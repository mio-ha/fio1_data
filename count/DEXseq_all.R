source("/home/miyokawa/ドキュメント/ngs/software/Subread_to_DEXSeq/load_SubreadOutput.R")
suppressPackageStartupMessages({
  require(dplyr)
})

# Load featureCounts result file
samp <- data.frame(row.names = c("Parker_Col-0_04c_1", "Parker_Col-0_04c_2", "Parker_Col-0_04c_3", "Parker_Col-0_04c_4", "Parker_Col-0_04c_5", "Parker_Col-0_04c_6", "Parker_Col-0_12c_1", "Parker_Col-0_12c_2", "Parker_Col-0_12c_3", "Parker_Col-0_12c_4", "Parker_Col-0_12c_5", "Parker_Col-0_12c_6", "Parker_Col-0_20c_1", "Parker_Col-0_20c_2", "Parker_Col-0_20c_3", "Parker_Col-0_20c_4", "Parker_Col-0_20c_5", "Parker_Col-0_20c_6", "Parker_Col-0_28c_1", "Parker_Col-0_28c_2", "Parker_Col-0_28c_3", "Parker_Col-0_28c_4", "Parker_Col-0_28c_5", "Parker_Col-0_28c_6", "Parker_fio1_04c_1", "Parker_fio1_04c_2", "Parker_fio1_04c_3", "Parker_fio1_04c_4", "Parker_fio1_04c_5", "Parker_fio1_04c_6", "Parker_fio1_12c_1", "Parker_fio1_12c_2", "Parker_fio1_12c_3", "Parker_fio1_12c_4", "Parker_fio1_12c_5", "Parker_fio1_12c_6", "Parker_fio1_20c_1", "Parker_fio1_20c_2", "Parker_fio1_20c_3", "Parker_fio1_20c_4", "Parker_fio1_20c_5", "Parker_fio1_20c_6", "Parker_fio1_28c_1", "Parker_fio1_28c_2", "Parker_fio1_28c_3", "Parker_fio1_28c_4", "Parker_fio1_28c_5", "Parker_fio1_28c_6"), 
                   condition = c("Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1"), 
                   type = c("04c", "04c", "04c", "04c", "04c", "04c", "12c", "12c", "12c", "12c", "12c", "12c", "20c", "20c", "20c", "20c", "20c", "20c", "28c", "28c", "28c", "28c", "28c", "28c", "04c", "04c", "04c", "04c", "04c", "04c", "12c", "12c", "12c", "12c", "12c", "12c", "20c", "20c", "20c", "20c", "20c", "20c", "28c", "28c", "28c", "28c", "28c", "28c"))
dxd.fc <- DEXSeqDataSetFromFeatureCounts("/media/miyokawa/8TB-Data3/tech_data/count/dex_counts_Parker.txt",
                                         flattenedfile = "Araport11_GTF_genes_transposons.Apr2023_flat.gtf",
                                         sampleData = samp)

# Adjust for coverage biases among samples
dxd <- estimateSizeFactors(dxd.fc)

# Estimate variance or dispersion parameters individually exon by exon
dxd <- estimateDispersions(dxd)

png(filename = "test.png", width = 600, height = 600)
plotDispEsts(dxd)
dev.off()

# Test for differential exon usage in each gene
dxd <- testForDEU(dxd)

# Estimate relative exon usage fold changes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")

# Output DEXSeq result
dxr1 <- DEXSeqResults(dxd)
table <- as.data.frame(dxr1)
table$transcripts <- as.character(table$transcripts)
write.table(table, quote = F, sep = "\t", file = "DEXSeq_Wang_result.txt", row.names=F)
