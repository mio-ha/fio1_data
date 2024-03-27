import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

study = ["Cai", "Parker20c", "Sun", "Wang"]

agg = pd.read_table("ag_FLOR_rMATS4.tsv", header=None, usecols=range(10))
merge = agg#[(agg["Score"] < 0.05) & (agg["Freq"] > 1)]
merge.columns = ["Chr", "Start", "End", "Geneid", "Start2", "End2", "logFC", "Strand", "Symbol", "Pathway"]
merge.drop("Pathway", axis=1, inplace=True)
merge["Name"] = merge["Geneid"] + ":" + merge["Start"].astype(str) + ":" + merge["End"].astype(str)

#index = len(study)
for i in study:
	sp = pd.read_table("rmats/ALL.MATS.JCEC_fix_FC_{}2.csv".format(i))
	#sp.columns = ["Chr", "Start", "End", "Geneid", "FDR", "Strand", "Diff", "abs"]
	sp.columns = ["Chr", "Start", "End", "Geneid", "Strand", "Diff", "abs"]
	sp["Name"] = sp["Geneid"] + ":" + sp["Start"].astype(str) + ":" + sp["End"].astype(str)

	sp2 = sp.loc[:, ["Name", "Diff"]]
	sp2.rename(columns={"Diff": i}, inplace=True)

	merge = pd.merge(merge, sp2, on="Name", how="left")
	
	p = sp.loc[:, ["Name", "Chr", "Start", "End", "Strand"]]
	if i == study[0]:
		pos = p
	else:
		pos = pd.concat([pos, p])

data = pd.concat([merge["Start"], merge.iloc[:, 8:14]], axis=1)
#data = merge.iloc[:, 8:14]
data["average"] = data.iloc[:, 3:7].mean(axis=1)
#data.iloc[:, 3:8] = -data.iloc[:, 3:8] #fold changeの正負が反転している場合

#data.set_index(merge["Name"], inplace=True)
data.sort_values("average", ascending=False, inplace=True)

pos.drop_duplicates(inplace=True)
out = pd.merge(pos, data, on="Name", how='right')
#out.to_csv("ag_logFC_rMATS3.txt", sep="\t", index=False)

data['count'] = data.groupby('Symbol').cumcount() + 1
name_counts = data['Symbol'].value_counts()
data['new_symbol'] = "$" + data['Symbol'] + r"\:_{" + data["Start"].astype(str) + "}$"
#data['new_symbol'] = data.apply(lambda row: row['Symbol'] if name_counts[row['Symbol']] == 1 else f"{row['count']}_{row['Symbol']}", axis=1)
max_label_len = max(len(label) for label in data["new_symbol"])-4
data.set_index(data["new_symbol"], inplace=True)
data.drop(["Start", "Symbol", "Name", "average", "count", "new_symbol"], axis=1, inplace=True)
data.columns = ["A", "B", "C", "D"]
print(data)

sns.set(font_scale=1.5)
#plt.rcParams["font.size"] = 12
sns.set_style(style='white')
#plt.rcParams['figure.subplot.bottom'] = 0.35
fig, ax = plt.subplots(figsize=(5, 9))
p = sns.heatmap(data, vmax=1, vmin=-1, cmap='coolwarm', square=True, yticklabels=False, linecolor='dimgrey', linewidth=1)
plt.tick_params(left=False)
plt.ylabel("")
colorbar = p.collections[0].colorbar
colorbar.set_label('Difference in PSI')
plt.subplots_adjust(left=0.3, right=0.88)

for i in range(len(data)):
    position = -0.15 * max_label_len
    ax.text(position, i+0.5, data.index[i], ha='left', va='center')

plt.savefig("meta_rMATS9.pdf")
#plt.show()
