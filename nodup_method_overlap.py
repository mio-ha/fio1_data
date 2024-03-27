import pandas as pd
import os

study = ["Cai", "Parker20c", "Sun", "Wang"]
method = ["SUPPA", "DEX", "rMATS"]

for i in study:
	for j in range(len(method)):
		for k in range(len(method)-1):
			if j >= k+1:
				continue

			input = "overlap_{}_{}_{}.bed".format(i, method[j], method[k+1])
			name = os.path.splitext(os.path.basename(input))[0]
			peak = pd.read_table(input, header=None)

			peak.drop_duplicates(inplace=True, subset=1)
			peak.to_csv("{}_nodup.bed".format(name), sep="\t", index=False, header=False)

