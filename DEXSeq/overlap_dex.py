import pandas as pd
import re
import os

file_list = ["DEXSeq_Cai_result.txt", "DEXSeq_Parker20c_result.txt", "DEXSeq_Sun2_result.txt", "DEXSeq_Wang_result.txt"]
list_n = len(file_list)

d = {}
for i in range(list_n):
	d[i] = pd.read_table(file_list[i], low_memory=False)

sig_d = {}
for i, df in d.items():
	df["uniqueID"] = df["groupID"] + ":" + df["featureID"]
	sig_df = df[df["padj"] < 0.05]
	sig_d[i] = sig_df

for i in range(list_n-1):
	for j in range(list_n-1):
		if i >= j+1:
			continue
			
		merge = pd.merge(sig_d[i], sig_d[j+1], on="uniqueID", how='inner')
		a = re.sub(".*_(.*?)_.*", "\\1", file_list[i])
		b = re.sub(".*_(.*?)_.*", "\\1", file_list[j+1])
		print(a, b, len(merge), "({}, {}: {})".format(len(sig_d[i]), len(sig_d[j+1]), len(merge)/min(len(sig_d[i]), len(sig_d[j+1]))))
		
		outfile = "overlap/DEXSeq_{}_{}_overlap.txt".format(a, b)
		is_file = os.path.isfile(outfile)
		if is_file:
			print("{} already exists.".format(outfile))
			continue
		merge.to_csv(outfile, sep="\t", index=False)
