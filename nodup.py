import pandas as pd
import os

input = "overlap/overlap_Wang_SUPPA_edgeR.bed"
name = os.path.splitext(os.path.basename(input))[0]
peak = pd.read_table(input, header=None)

peak.drop_duplicates(inplace=True)
peak.to_csv("overlap/{}_nodup.bed".format(name), sep="\t", index=False, header=False)

