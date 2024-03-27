import pandas as pd
import os
import sys

study = ["Cai", "Parker20c", "Sun", "Wang"]

list_n = len(study)

for i in range(list_n-1):
	for j in range(list_n-1):
		if i >= j+1:
			continue

		outfile = "overlap_SUPPA_{}_{}.txt".format(study[i], study[j+1])
		is_file = os.path.isfile(outfile)
		if is_file:
			print("{} already exists.".format(outfile))
			continue

		df = pd.read_table("SUPPA/ara_diffSplice_{}.dpsi".format(study[i]), sep = "[;\t]", engine='python', header=None, skiprows=1)
		df2 = pd.read_table("SUPPA/ara_diffSplice_{}.dpsi".format(study[j+1]), sep = "[;\t]", engine='python', header=None, skiprows=1)
		df.columns = ["gene", "event", "dpsi_{}".format(study[i]), "p_val_{}".format(study[i])]
		df2.columns = ["gene", "event", "dpsi_{}".format(study[j+1]), "p_val_{}".format(study[j+1])]
		sig_df = df[df["p_val_{}".format(study[i])] < 0.05]
		sig_df2 = df2[df2["p_val_{}".format(study[j+1])] < 0.05]
		print(study[i], len(sig_df), study[j+1], len(sig_df2))
		#print(sig_df.columns)
		overlap = pd.merge(sig_df, sig_df2, on=["gene", "event"], how = 'inner')
		print(study[i], study[j+1], len(overlap))
		#sys.exit()
		overlap.to_csv(outfile, sep="\t", index=False)
			
