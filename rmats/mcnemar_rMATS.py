import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar

file_list_a = np.array(["rmats_Parker20c/ALL.MATS.JCEC.csv", "rmats_Parker_100bp/ALL.MATS.JCEC.csv", "rmats_Parker_50bp/ALL.MATS.JCEC.csv"])
file_list_a2 = np.roll(file_list_a, 1)

file_list_b = np.char.replace(file_list_a, "MATS.JCEC", "MATS.NS.JCEC")
file_list_b2 = np.roll(file_list_b, 1)

for i in range(len(file_list_a)):
	a = pd.read_table(file_list_a[i], low_memory=False)
	b = pd.read_table(file_list_a2[i], low_memory=False)
	c = pd.read_table(file_list_b[i], low_memory=False)
	d = pd.read_table(file_list_b2[i], low_memory=False)
	print(file_list_a[i], file_list_a2[i])

	a_S = a
	a_NS = c
	b_S = b
	b_NS = d

	aSbS = pd.merge(a_S, b_S, on=["GeneID", "longExonStart_0base"], how='inner')
	aSbNS = pd.merge(a_S, b_NS, on=["GeneID", "longExonStart_0base"], how='inner')
	aNSbS = pd.merge(a_NS, b_S, on=["GeneID", "longExonStart_0base"], how='inner')
	aNSbNS = pd.merge(a_NS, b_NS, on=["GeneID", "longExonStart_0base"], how='inner')

	obs = [[len(aSbS), len(aSbNS)], [len(aNSbS), len(aNSbNS)]]
	print(obs)

	result = mcnemar(obs, exact=True)
	print("McNemar's test statistic:", result.statistic)
	print("p-value:", result.pvalue)
