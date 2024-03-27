import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar

file_list = np.array(["ara_diffSplice_Parker20c.dpsi", "ara_diffSplice_Parker20c_4rep.dpsi", "ara_diffSplice_Parker20c_2rep.dpsi"])
file_list2 = np.roll(file_list, 1)

for i in range(len(file_list)):
	a = pd.read_table(file_list[i], low_memory=False)
	b = pd.read_table(file_list2[i], low_memory=False)
	a.reset_index(inplace=True)
	b.reset_index(inplace=True)
	print(file_list[i], file_list2[i])

	a_S = a[a.iloc[:, 2] < 0.05]
	a_NS = a[a.iloc[:, 2] >= 0.05]
	b_S = b[b.iloc[:, 2] < 0.05]
	b_NS = b[b.iloc[:, 2] >= 0.05]

	aSbS = pd.merge(a_S, b_S, on="index", how='inner')
	aSbNS = pd.merge(a_S, b_NS, on="index", how='inner')
	aNSbS = pd.merge(a_NS, b_S, on="index", how='inner')
	aNSbNS = pd.merge(a_NS, b_NS, on="index", how='inner')

	obs = [[len(aSbS), len(aSbNS)], [len(aNSbS), len(aNSbNS)]]
	print(obs)

	result = mcnemar(obs, exact=True)
	print("McNemar's test statistic:", result.statistic)
	print("p-value:", result.pvalue)
