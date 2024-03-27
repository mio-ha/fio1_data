import pandas as pd
import os
import sys

dir = "rmats_Cai"
event_list = ["A3SS", "A5SS", "MXE", "RI", "SE"]
EC = True

if EC == True:
	saffix = "JCEC"
else:
	saffix = "JC"
outfile = "{}/ALL.MATS.{}.txt".format(dir, saffix)
is_file = os.path.isfile(outfile)
if is_file:
	print("{} already exists.".format(outfile))
	sys.exit()

total = 0
print(dir)

for ev in event_list:
	df = pd.read_table("{}/{}.MATS.{}.txt".format(dir, ev, saffix))
	sig_df = df[df["FDR"] < 0.05]
	total += len(sig_df)
	print(ev, len(sig_df))
	sig_df.to_csv(outfile, mode='a', sep="\t", index=False)

print("total", total)
