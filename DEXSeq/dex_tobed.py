import pandas as pd

openfile = "DEXSeq_Wang_result.txt"
outfile = "DEXSeq_Wang.bed"

df = pd.read_table(openfile, low_memory=False)
padj_df = df[df["padj"] < 0.05]
bed = pd.concat([padj_df["genomicData.seqnames"], padj_df["genomicData.start"], padj_df["genomicData.end"], padj_df["groupID"], padj_df["padj"], padj_df["genomicData.strand"]], axis=1)
bed.to_csv(outfile, sep="\t", index=False, header=False)
