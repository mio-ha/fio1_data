import pandas as pd
import os

study = ["DEXSeq", "edgeR", "rMATS", "SUPPA"]
method = "Cai"
n_study = len(study)

for i in range(n_study):
	for j in range(n_study):
		#if i >= j+1:
		#	continue

		input = "overlap/overlap_reverse/overlap_fix_{}_{}_{}.bed".format(method, study[i], study[j])
		if os.path.isfile(input) == False:
			continue
		name = os.path.splitext(os.path.basename(input))[0]
		peak = pd.read_table(input, header=None)

		peak.drop_duplicates(inplace=True)
		peak.to_csv("overlap/overlap_reverse/{}_nodup.bed".format(name), sep="\t", index=False, header=False)

