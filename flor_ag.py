import pandas as pd

#flor = pd.read_table("/home/miyokawa/ダウンロード/FLOR-ID.txt", header=None, sep="\n")
with open("/home/miyokawa/ダウンロード/FLOR-ID.txt", 'r') as file:
	flor = file.readlines()
	flor = [line.strip() for line in flor]
ag_df = pd.read_table("ag_logFC_rMATS3.txt", header=None, skiprows=1)

flo_df = ag_df[ag_df[0].str.contains('|'.join(flor))]
print(flo_df)
flo_df.to_csv("ag_FLOR_rMATS4.tsv", sep="\t", header=False, index=False)
