import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar

a = pd.read_table("DEXSeq_Parker20c_2rep_result.txt", low_memory=False)
b = pd.read_table("DEXSeq_Parker20c_4rep_result.txt", low_memory=False)

a_S = a[a["padj"] < 0.05]
a_NS = a[a["padj"] >= 0.05]
b_S = b[b["padj"] < 0.05]
b_NS = b[b["padj"] >= 0.05]

aSbS = pd.merge(a_S, b_S, on=["groupID", "featureID"], how='inner')
aSbNS = pd.merge(a_S, b_NS, on=["groupID", "featureID"], how='inner')
aNSbS = pd.merge(a_NS, b_S, on=["groupID", "featureID"], how='inner')
aNSbNS = pd.merge(a_NS, b_NS, on=["groupID", "featureID"], how='inner')

obs = [[len(aSbS), len(aSbNS)], [len(aNSbS), len(aNSbNS)]]
print(obs)

result = mcnemar(obs, exact=True)
print("McNemar's test statistic:", result.statistic)
print("p-value:", result.pvalue)
