source("/home/miyokawa/ドキュメント/ngs/software/Subread_to_DEXSeq/load_SubreadOutput.R")
suppressPackageStartupMessages({
  require(dplyr)
})

rn <- c("Cai_Col-0_1.bam", "Cai_Col-0_2.bam", "Cai_Col-0_3.bam", "Cai_fio1_1.bam", "Cai_fio1_2.bam", "Cai_fio1_3.bam", "Parker_Col-0_20c_1.bam", "Parker_Col-0_20c_2.bam", "Parker_Col-0_20c_3.bam", "Parker_Col-0_20c_4.bam", "Parker_Col-0_20c_5.bam", "Parker_Col-0_20c_6.bam", "Parker_fio1_20c_1.bam", "Parker_fio1_20c_2.bam", "Parker_fio1_20c_3.bam", "Parker_fio1_20c_4.bam", "Parker_fio1_20c_5.bam", "Parker_fio1_20c_6.bam", "Sun_Col-0_1.bam", "Sun_Col-0_2.bam", "Sun_fio1-1_1.bam", "Sun_fio1-1_2.bam", "Sun_fio1-5_1.bam", "Sun_fio1-5_2.bam", "Wang_Col-0_1.bam", "Wang_Col-0_2.bam", "Wang_fio1_1.bam", "Wang_fio1_2.bam")
condition <- factor(c("Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1", "fio1", "fio1", "Col-0", "Col-0", "fio1", "fio1"))
type <- factor(c("Cai", "Cai", "Cai", "Cai", "Cai", "Cai", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker", "Parker","Parker", "Sun", "Sun", "Sun", "Sun", "Sun", "Sun", "Wang", "Wang", "Wang", "Wang"))
samp <- data.frame(row.names=rn, condition, type)
dxd.fc <- DEXSeqDataSetFromFeatureCounts("/media/miyokawa/8TB-Data3/tech_data/count/dex_counts_fio1_20c.txt",
                                         flattenedfile = "../count/Araport11_GTF_genes_transposons.Apr2023_flat.gtf",
                                         sampleData = samp, design= ~ sample + exon + condition:exon)

dxd <- estimateSizeFactors(dxd.fc)

formulaFullModel    =  ~ sample + exon + type:exon + condition:exon
formulaReducedModel =  ~ sample + exon + type:exon 

multicoreParam <- MulticoreParam(workers = 8)
dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=multicoreParam)

# Test for differential exon usage in each gene
dxd = testForDEU( dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel, BPPARAM=multicoreParam)

# Estimate relative exon usage fold changes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=multicoreParam)

dxr1 <- DEXSeqResults(dxd)
table <- as.data.frame(dxr1)
table$transcripts <- as.character(table$transcripts)
write.table(table, quote = F, sep = "\t", file = "DEXSeq_all_rem_batch_result.txt", row.names=F)
