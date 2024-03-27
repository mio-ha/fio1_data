source("/home/miyokawa/ドキュメント/ngs/software/Subread_to_DEXSeq/load_SubreadOutput.R")
suppressPackageStartupMessages({
  require(dplyr)
})

rn <- c("Cai_Col-0_1.bam", "Cai_Col-0_2.bam", "Cai_Col-0_3.bam", "Cai_fio1_1.bam", "Cai_fio1_2.bam", "Cai_fio1_3.bam", "Parker_Col-0_20c_1.bam", "Parker_Col-0_20c_2.bam", "Parker_Col-0_20c_3.bam", "Parker_Col-0_20c_4.bam", "Parker_Col-0_20c_5.bam", "Parker_Col-0_20c_6.bam", "Parker_fio1_20c_1.bam", "Parker_fio1_20c_2.bam", "Parker_fio1_20c_3.bam", "Parker_fio1_20c_4.bam", "Parker_fio1_20c_5.bam", "Parker_fio1_20c_6.bam", "Sun_Col-0_1.bam", "Sun_Col-0_2.bam", "Sun_fio1-1_1.bam", "Sun_fio1-1_2.bam", "Sun_fio1-5_1.bam", "Sun_fio1-5_2.bam", "Wang_Col-0_1.bam", "Wang_Col-0_2.bam", "Wang_fio1_1.bam", "Wang_fio1_2.bam")
condition <- factor(c("Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1"))
samp <- data.frame(row.names=rn, condition)
dxd.fc <- DEXSeqDataSetFromFeatureCounts("/media/miyokawa/8TB-Data3/tech_data/count/dex_counts_fio1_20c.txt",
                                         flattenedfile = "../count/Araport11_GTF_genes_transposons.Apr2023_flat.gtf",
                                         sampleData = samp)

dxd <- estimateSizeFactors(dxd.fc)

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
write.table(table, quote = F, sep = "\t", file = "../DEXSeq/DEXSeq_all_norem_batch_result.txt", row.names = F)
