import numpy as np
import pandas as pd

type = ["A3SS", "A5SS", "MXE", "RI", "SE"]

for i in type:
	adj = pd.read_table("rmats_all/{}_combat_adj.txt".format(i), header=None)
	adj = adj.astype(np.int64)
	ijc1 = adj.iloc[:, 1:14]
	print(ijc1)
