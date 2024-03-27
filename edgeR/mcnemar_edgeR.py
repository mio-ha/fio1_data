import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar

file_list = np.array(["edgeR_result_splicing_Parker.txt", "edgeR_result_splicing_Parker_100bp.txt", "edgeR_result_splicing_Parker_50bp.txt"])
file_list2 = np.roll(file_list, 1)

for i in range(len(file_list)):
	a = pd.read_table(file_list[i], low_memory=False)
	b = pd.read_table(file_list2[i], low_memory=False)
	print(file_list[i], file_list2[i])

	a_S = a[a["FDR"] < 0.05]
	a_NS = a[a["FDR"] >= 0.05]
	b_S = b[b["FDR"] < 0.05]
	b_NS = b[b["FDR"] >= 0.05]

	aSbS = pd.merge(a_S, b_S, on=["Geneid", "Start"], how='inner')
	aSbNS = pd.merge(a_S, b_NS, on=["Geneid", "Start"], how='inner')
	aNSbS = pd.merge(a_NS, b_S, on=["Geneid", "Start"], how='inner')
	aNSbNS = pd.merge(a_NS, b_NS, on=["Geneid", "Start"], how='inner')

	obs = [[len(aSbS), len(aSbNS)], [len(aNSbS), len(aNSbNS)]]
	print(obs)

	result = mcnemar(obs, exact=True)
	print("McNemar's test statistic:", result.statistic)
	print("p-value:", result.pvalue)
