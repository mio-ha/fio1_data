import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

study = ["Cai", "Parker", "Sun", "Wang"]

agg = pd.read_csv("ag_score_DEXSeq3.txt", sep='\t|:', engine="python", header=None, skiprows=1)
agg.columns = ["Name", "Start", "End", "Score", "Freq"]
merge = agg[(agg["Score"] < 0.05) & (agg["Freq"] > 1)]

#index = len(study)
for i in study:
	sp = pd.read_table("DEXSeq/DEXSeq_{}_sig_abs.csv".format(i))
	#sp["Name"] = sp["groupID"] + ":" + sp["featureID"]
	sp2 = sp.loc[:, ["groupID", "genomicData.start", "genomicData.end", "log2fold_fio1_Col_0"]]
	#sp2.rename(columns={"log2fold_fio1_Col_0": i}, inplace=True)
	sp2.columns = ["Name", "Start", "End", i]
	merge = pd.merge(merge, sp2, on=["Name", "Start", "End"], how="left")
	
	p = sp.loc[:, ["groupID", "featureID", "genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand"]]
	if i == study[0]:
		pos = p
	else:
		pos = pd.concat([pos, p])

print(merge)
data = merge.iloc[:, 5:9]
data["average"] = data.mean(axis=1)
data = pd.concat([merge["Name"], merge["Start"], merge["End"], merge["Score"], data], axis=1)
data.sort_values("average", ascending=False, inplace=True)

pos.rename(columns={"groupID": "Name", "genomicData.seqnames": "Chr", "genomicData.start": "Start", "genomicData.end": "End", "genomicData.strand": "Strand"}, inplace=True)
pos.drop_duplicates(inplace=True)
out = pd.merge(pos, data, on=["Name", "Start", "End"], how='right')
out.drop("average", axis=1, inplace=True)
print(out)
out.to_csv("ag_logFC_DEXSeq_all.txt", sep="\t", index=False)
sys.exit()

data.drop("average", axis=1, inplace=True)
sns.set(font_scale=2)

#plt.rcParams["font.size"] = 12
sns.set_style(style='white')
#plt.rcParams['figure.subplot.bottom'] = 0.35
fig, ax = plt.subplots(figsize=(9, 9))
p = sns.heatmap(data, vmax=5, vmin=-5, cmap='coolwarm')
plt.tick_params(labelleft=False, left=False)
plt.ylabel("")

plt.savefig("meta_DEXSeq2.pdf")

