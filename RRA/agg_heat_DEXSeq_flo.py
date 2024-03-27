import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

study = ["Cai", "Parker", "Sun", "Wang"]

agg = pd.read_table("ag_FLOR_DEXSeq2.tsv", header=None)
merge = agg#[(agg["Score"] < 0.05) & (agg["Freq"] > 1)]
merge.columns = ["Chr", "Start", "End", "Geneid", "Exon", "logFC", "Strand", "Symbol", "Pathway"]
merge.drop("Pathway", axis=1, inplace=True)
merge["Name"] = merge["Geneid"] + ":" + merge["Exon"]

#index = len(study)
for i in study:
	sp = pd.read_table("DEXSeq/DEXSeq_{}_sig_abs.csv".format(i))
	sp["Name"] = sp["groupID"] + ":" + sp["featureID"]
	sp2 = sp.loc[:, ["Name", "log2fold_fio1_Col_0"]]
	sp2.rename(columns={"log2fold_fio1_Col_0": i}, inplace=True)
	merge = pd.merge(merge, sp2, on="Name", how="left")
	
	p = sp.loc[:, ["Name", "genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand"]]
	if i == study[0]:
		pos = p
	else:
		pos = pd.concat([pos, p])

#merge["Symbol"] = merge["Start"].astype(str) +": "+ merge["Symbol"]
merge["Symbol"] = "$" + merge['Symbol'] + r"\:_{" + merge["Start"].astype(str) + "}$"
data = merge.iloc[:, 7:13]
data["average"] = data.iloc[:, 2:6].mean(axis=1)
#data.set_index(merge["Name"], inplace=True)
data.sort_values("average", ascending=False, inplace=True)

pos.rename(columns={"genomicData.seqnames": "chr", "genomicData.start": "start", "genomicData.end": "end", "genomicData.strand": "strand"}, inplace=True)
pos.drop_duplicates(inplace=True)
out = pd.merge(pos, data, on="Name", how='right')
out.to_csv("ag_logFC_DEXSeq5.txt", sep="\t", index=False)
sys.exit()

data['count'] = data.groupby('Symbol').cumcount() + 1
name_counts = data['Symbol'].value_counts()
data['new_symbol'] = data["Symbol"]#data.apply(lambda row: row['Symbol'] if name_counts[row['Symbol']] == 1 else f"{row['count']}_{row['Symbol']}", axis=1)
max_label_len = max(len(label) for label in data["new_symbol"])-4
data.set_index(data["new_symbol"], inplace=True)
data.drop(["Symbol", "Name", "average", "count", "new_symbol"], axis=1, inplace=True)
data = data.iloc[:-1] #FIO1を削除
data.columns = ["A", "B", "C", "D"]

sns.set(font_scale=1.1)
#plt.rcParams["font.size"] = 12
sns.set_style(style='white')
#plt.rcParams['figure.subplot.bottom'] = 0.35
fig, ax = plt.subplots(figsize=(5, 9))
p = sns.heatmap(data, vmax=5, vmin=-5, cmap='coolwarm', square=True, yticklabels=False, linecolor='dimgrey', linewidth=1)
plt.tick_params(left=False)
plt.ylabel("")
colorbar = p.collections[0].colorbar
colorbar.set_label('log2(FC)')

for i in range(len(data)):
    position = -0.36 * max_label_len
    ax.text(position, i+0.5, data.index[i], ha='left', va='center')

plt.savefig("meta_DEXSeq6.pdf", dpi=300)
#plt.show()
