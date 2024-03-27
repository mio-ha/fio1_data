import pandas as pd
import os

study = ["Cai", "Parker", "Sun", "Wang"]
method = "SUPPA"
n_study = len(study)

for i in range(n_study-1):
	for j in range(n_study-1):
		if i >= j+1:
			continue

		input = "overlap/{}_{}_{}_overlap2.bed".format(method, study[i], study[j+1])
		name = os.path.splitext(os.path.basename(input))[0]
		peak = pd.read_table(input, header=None)

		peak.drop_duplicates(inplace=True)
		peak.to_csv("overlap/{}_nodup.bed".format(name), sep="\t", index=False, header=False)

