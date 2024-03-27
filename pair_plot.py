import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

cai = pd.read_table("DEXSeq/DEXSeq_Cai_sig_abs.csv", usecols=[0, 1, 9])
parker = pd.read_table("DEXSeq/DEXSeq_Parker_sig_abs.csv", usecols=[0, 1, 9])
sun = pd.read_table("DEXSeq/DEXSeq_Sun_sig_abs.csv", usecols=[0, 1, 9])
wang = pd.read_table("DEXSeq/DEXSeq_Wang_sig_abs.csv", usecols=[0, 1, 9])

cai.rename(columns={"log2fold_fio1_Col_0": "Cai"}, inplace=True)
parker.rename(columns={"log2fold_fio1_Col_0": "Parker"}, inplace=True)
sun.rename(columns={"log2fold_fio1_Col_0": "Sun"}, inplace=True)
wang.rename(columns={"log2fold_fio1_Col_0": "Wang"}, inplace=True)

merge = pd.merge(cai, parker, on=["groupID", "featureID"], how="outer")
merge = pd.merge(merge, sun, on=["groupID", "featureID"], how="outer")
merge = pd.merge(merge, wang, on=["groupID", "featureID"], how="outer")

merge.to_csv("log2_4study.csv", sep="\t")
pg = sns.pairplot(merge, kind='reg', plot_kws={"scatter_kws": {"s": 3, "alpha": 0.5}, 'ci': None,  "truncate": False, "line_kws": {"color": "grey", "linewidth": 1}})

for i in range(4):
	for j in range(4):
		if i == j:
			continue
		model = LinearRegression()
		nona = merge.dropna(subset=[merge.columns[2+i], merge.columns[2+j]])
		x = np.array(nona.iloc[:, 2+i]).reshape(-1,1)
		y = np.array(nona.iloc[:, 2+j]).reshape(-1,1)
		model.fit(x, y)
		y_pred = model.predict(x)
		r2 = r2_score(y, y_pred)
		print(merge.columns[2+i], merge.columns[2+j], r2)

#print(merge.iloc[:, 2:])
#hm = sns.heatmap(merge.iloc[:, 2:].corr(), annot=True)

pg.axes[0,0].set_xlim((-5,5))
pg.axes[0,1].set_xlim((-5,5))
pg.axes[0,2].set_xlim((-5,5))
pg.axes[0,3].set_xlim((-5,5))
pg.axes[0,0].set_ylim((-5,5))
pg.axes[1,0].set_ylim((-5,5))
pg.axes[2,0].set_ylim((-5,5))
pg.axes[3,0].set_ylim((-5,5))

plt.show()
