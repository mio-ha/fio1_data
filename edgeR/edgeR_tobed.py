import pandas as pd

read_file = "edgeR_result_splicing_Cai.txt"
out_file = "edgeR_result_splicing_Cai.bed"

df = pd.read_table(read_file, low_memory=False)
sig = df[df["FDR"] < 0.05]
bed = sig.reindex(columns=["Chr", "Start", "End", "Geneid", "FDR", "Strand"])
bed.to_csv(out_file, sep="\t", index=False, header=False)
