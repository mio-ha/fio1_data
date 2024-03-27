library(RobustRankAggreg)
set.seed(1)

mode <- 4 #DEXSeq:1, edgeR:2, rMATS:3, SUPPA4

#DEXSeq
if (mode == 1) {
	cai = read.table("DEXSeq/DEXSeq_Cai_sig_abs.csv", sep="\t", header=T)
	parker = read.table("DEXSeq/DEXSeq_Parker_sig_abs.csv", sep="\t", header=T)
	sun = read.table("DEXSeq/DEXSeq_Sun_sig_abs.csv", sep="\t", header=T)
	wang = read.table("DEXSeq/DEXSeq_Wang_sig_abs.csv", sep="\t", header=T)

	cai_all = read.table("DEXSeq/DEXSeq_Cai_result.txt", sep="\t", header=T)
	parker_all = read.table("DEXSeq/DEXSeq_Parker20c_result.txt", sep="\t", header=T)
	sun_all = read.table("DEXSeq/DEXSeq_Sun2_result.txt", sep="\t", header=T)
	wang_all = read.table("DEXSeq/DEXSeq_Wang_result.txt", sep="\t", header=T)
}

#edgeR
if (mode == 2) {
	cai = read.table("edgeR/edgeR_result_splicing_Cai_sig.csv", sep="\t", header=T)
	parker = read.table("edgeR/edgeR_result_splicing_Parker_sig.csv", sep="\t", header=T)
	sun = read.table("edgeR/edgeR_result_splicing_Sun_sig.csv", sep="\t", header=T)
	wang = read.table("edgeR/edgeR_result_splicing_Wang_sig.csv", sep="\t", header=T)

	cai_all = read.table("edgeR/edgeR_result_splicing_Cai.txt", sep="\t", header=T)
	parker_all = read.table("edgeR/edgeR_result_splicing_Parker.txt", sep="\t", header=T)
	sun_all = read.table("edgeR/edgeR_result_splicing_Sun.txt", sep="\t", header=T)
	wang_all = read.table("edgeR/edgeR_result_splicing_Wang.txt", sep="\t", header=T)
}

#rMATS
if (mode == 3) {
	cai = read.table("rmats/ALL.MATS.JCEC_fix_FC_Cai2.csv", sep="\t", header=F)
	parker = read.table("rmats/ALL.MATS.JCEC_fix_FC_Parker20c2.csv", sep="\t", header=F)
	sun = read.table("rmats/ALL.MATS.JCEC_fix_FC_Sun2.csv", sep="\t", header=F)
	wang = read.table("rmats/ALL.MATS.JCEC_fix_FC_Wang2.csv", sep="\t", header=F)
}

#SUPPA
if (mode == 4) {
	cai = read.table("SUPPA/diffSplice/ara_diffSplice_Cai_fix_FC_uniq.csv", sep="\t", header=F)
	parker = read.table("SUPPA/diffSplice/ara_diffSplice_Parker20c_fix_FC_uniq.csv", sep="\t", header=F)
	sun = read.table("SUPPA/diffSplice/ara_diffSplice_Sun_fix_FC_uniq.csv", sep="\t", header=F)
	wang = read.table("SUPPA/diffSplice/ara_diffSplice_Wang_fix_FC_uniq.csv", sep="\t", header=F)
	
	cai_all = read.table("SUPPA/psi/Cai_event_psi_col0.txt", sep="\t", header=T)
	parker_all = read.table("SUPPA/psi/Parker20c_event_psi_col0.txt", sep="\t", header=T)
	sun_all = read.table("SUPPA/psi/Sun_event_psi_col0.txt", sep="\t", header=T)
	wang_all = read.table("SUPPA/psi/Wang_event_psi_col0.txt", sep="\t", header=T)
}

#DEXSeq
if (mode == 1) {
	cai_id <- paste(cai$groupID, cai$genomicData.start, sep=":")
	parker_id <- paste(parker$groupID, parker$genomicData.start, sep=":")
	sun_id <- paste(sun$groupID, sun$genomicData.start, sep=":")
	wang_id <- paste(wang$groupID, wang$genomicData.start, sep=":")
	cai_id <- paste(cai_id, cai$genomicData.end, sep=":")
	parker_id <- paste(parker_id, parker$genomicData.end, sep=":")
	sun_id <- paste(sun_id, sun$genomicData.end, sep=":")
	wang_id <- paste(wang_id, wang$genomicData.end, sep=":")
}

#edgeR
if (mode == 2) {
	cai_id <- paste(cai$Geneid, cai$Start, sep=":")
	parker_id <- paste(parker$Geneid, parker$Start, sep=":")
	sun_id <- paste(sun$Geneid, sun$Start, sep=":")
	wang_id <- paste(wang$Geneid, wang$Start, sep=":")
	cai_id <- paste(cai_id, cai$End, sep=":")
	parker_id <- paste(parker_id, parker$End, sep=":")
	sun_id <- paste(sun_id, sun$End, sep=":")
	wang_id <- paste(wang_id, wang$End, sep=":")
}

#rMATS/SUPPA
if (mode == 3 | mode == 4) {
	cai_id <- paste(cai$V4, cai$V2, sep=":")
	parker_id <- paste(parker$V4, parker$V2, sep=":")
	sun_id <- paste(sun$V4, sun$V2, sep=":")
	wang_id <- paste(wang$V4, wang$V2, sep=":")
	cai_id <- paste(cai_id, cai$V3, sep=":")
	parker_id <- paste(parker_id, parker$V3, sep=":")
	sun_id <- paste(sun_id, sun$V3, sep=":")
	wang_id <- paste(wang_id, wang$V3, sep=":")
}

glist <- list(cai_id, parker_id, sun_id, wang_id)
if (mode == 1 | mode == 2) {
	num <- max(c(nrow(cai_all), nrow(parker_all), nrow(sun_all), nrow(wang_all)))
}
if (mode == 3) {
	num <- 15883
}
if (mode == 4) {
	num <- max(c(nrow(cai_all), nrow(parker_all), nrow(sun_all), nrow(wang_all)))
}
r = rankMatrix(glist, N=num)

freq=as.data.frame(table(unlist(glist)))
#ag=aggregateRanks(glist = glist)
ag = aggregateRanks(rmat = r)
ag$Freq=freq[match(ag$Name,freq$Var1),2]

if (mode == 1) {
	write.table(ag, "ag_score_DEXSeq4.txt", sep="\t", row.names=F)
}
if (mode == 2) {
	write.table(ag, "ag_score_edgeR4.txt", sep="\t", row.names=F)
}
if (mode == 3) {
	write.table(ag, "ag_score_rMATS4.txt", sep="\t", row.names=F)
}
if (mode == 4) {
	write.table(ag, "ag_score_SUPPA4.txt", sep="\t", row.names=F)
}
