import pandas as pd
import os
import sys

study = ["Cai", "Parker20c", "Sun", "Wang"]
event_list = ["A3SS", "A5SS", "MXE", "RI", "SE"]
start_list = ["longExonStart_0base", "longExonStart_0base", "1stExonStart_0base", "riExonStart_0base", "exonStart_0base"]
end_list = ["longExonEnd", "longExonEnd", "1stExonEnd", "riExonEnd", "exonEnd"]
EC = True

if EC == True:
	saffix = "JCEC"
else:
	saffix = "JC"

list_n = len(study)
nomerge = True

for i in range(list_n-1):
	for j in range(list_n-1):
		if i >= j+1:
			continue

		outfile = "overlap_rMATS.{}.{}_{}.txt".format(saffix, study[i], study[j+1])
		is_file = os.path.isfile(outfile)
		if is_file:
			print("{} already exists.".format(outfile))
			continue

		for ev_num in range(len(event_list)):
			df = pd.read_table("rmats_{}/{}.MATS.{}.txt".format(study[i], event_list[ev_num], saffix))
			df2 = pd.read_table("rmats_{}/{}.MATS.{}.txt".format(study[j+1], event_list[ev_num], saffix))
			sig_df = df[df["FDR"] < 0.05]
			sig_df2 = df2[df2["FDR"] < 0.05]
			#print(sig_df.columns)
			overlap = pd.merge(sig_df, sig_df2, on=["GeneID", "chr", "strand", start_list[ev_num], end_list[ev_num]], how = 'inner')
			#print(overlap)
			overlap.to_csv(outfile, mode='a', sep="\t", index=False)
			
