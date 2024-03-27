import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

study = ["Cai", "Parker", "Sun", "Wang"]

d = {}
nomerge = True

for i in range(len(study)):
	d[i] = {}
	for j in range(i+1, len(study)):
		dex = "overlap/DEXSeq_{}_{}_overlap.txt".format(study[i], study[j])
		try:
			d[i][j] = pd.read_table(dex, usecols=["log2fold_fio1_Col_0", "uniqueID"])
		except ValueError:
			pass
		d[i][j].rename(columns={"log2fold_fio1_Col_0": "log2fold_fio1_Col_0_{}_{}".format(study[i], study[j])})
		if nomerge == True:
			merge = d[i][j]
			nomerge = False
		else:
			merge = pd.merge(merge, d[i][j], on="uniqueID")

print(merge)
'''
df = pd.read_table("overlap/fc_table.csv")
corr = df.corr()
print(corr)

fig, ax = plt.subplots(figsize=(9, 9))
sns.heatmap(corr, vmax=0.5, cmap='Blues', annot=True, fmt="1.3f", square=True)
plt.show()
'''
