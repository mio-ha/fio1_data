import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

study = ["Cai", "Parker", "Sun", "Wang"]

agg = pd.read_table("ag_score_edgeR3.txt", sep='\t|:', engine="python", header=None, skiprows=1)
agg.columns = ["Name", "Start", "End", "Score", "Freq"]
merge = agg[(agg["Score"] < 0.05) & (agg["Freq"] > 1)]

#index = len(study)
for i in study:
	sp = pd.read_table("edgeR/edgeR_result_splicing_{}_sig.csv".format(i))
	sp.rename(columns={"Geneid": "Name"}, inplace=True)
	#sp["Name"] = sp["Geneid"] + ":" + sp["Start"].astype(str)
	sp2 = sp.loc[:, ["Name", "Start", "End", "logFC"]]
	sp2.rename(columns={"logFC": i}, inplace=True)
	print(sp2, merge)
	merge = pd.merge(merge, sp2, on=["Name", "Start", "End"], how="left")
	
	p = sp.loc[:, ["Name", "Chr", "Start", "End", "Strand"]]
	if i == study[0]:
		pos = p
	else:
		pos = pd.concat([pos, p])

data = merge.iloc[:, 5:9]
data["average"] = data.mean(axis=1)
data = pd.concat([merge["Name"], merge["Start"], merge["End"], merge["Score"], data], axis=1)
data.sort_values("average", ascending=False, inplace=True)

pos.drop_duplicates(inplace=True)
print(data, pos)
out = pd.merge(pos, data, on=["Name", "Start", "End"], how='right')
out.drop("average", axis=1, inplace=True)
print(out)
out.to_csv("ag_logFC_edgeR_all.txt", sep="\t", index=False)
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

plt.savefig("meta_edgeR2.pdf")

