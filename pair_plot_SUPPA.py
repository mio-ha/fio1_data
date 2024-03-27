import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

cai = pd.read_table("SUPPA/SUPPA_Cai_sig2.txt", usecols=[1, 3, 6], header=None)
parker = pd.read_table("SUPPA/SUPPA_Parker_sig2.txt", usecols=[1, 3, 6], header=None)
sun = pd.read_table("SUPPA/SUPPA_Sun_sig2.txt", usecols=[1, 3, 6], header=None)
wang = pd.read_table("SUPPA/SUPPA_Wang_sig2.txt", usecols=[1, 3, 6], header=None)
print(cai)
cai.rename(columns={6: "Cai"}, inplace=True)
parker.rename(columns={6: "Parker"}, inplace=True)
sun.rename(columns={6: "Sun"}, inplace=True)
wang.rename(columns={6: "Wang"}, inplace=True)

merge = pd.merge(cai, parker, on=[1, 3], how="outer")
merge = pd.merge(merge, sun, on=[1, 3], how="outer")
merge = pd.merge(merge, wang, on=[1, 3], how="outer")

merge.to_csv("log2_4study_rMATS.csv", sep="\t")


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

merge.drop(1, axis=1, inplace=True)
pg = sns.pairplot(merge, kind='reg', plot_kws={"scatter_kws": {"s": 3, "alpha": 0.5}, 'ci': None, "line_kws": {"color": "grey", "linewidth": 1}})

#print(merge.iloc[:, 1:])
#hm = sns.heatmap(merge.iloc[:, 1:].corr(), annot=True)


for i in range(4):
	for j in range(4):
		pg.axes[i, j].set_xlim((-1,1))
		pg.axes[i ,j].set_ylim((-1,1))

pg.map_offdiag(sns.regplot, truncate=False, line_kws={"color": "grey", "linewidth": 1}, scatter_kws={"s": 3, "alpha": 0.5}, ci=None)

plt.show()
