import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_table("SUPPA/overlap/overlap_table3.csv", index_col=0)
sns.set(font_scale=2)

#plt.rcParams["font.size"] = 12
sns.set_style(style='white')
#plt.rcParams['figure.subplot.bottom'] = 0.35
fig, ax = plt.subplots(figsize=(9, 9))
p = sns.heatmap(df, vmax=1, vmin=0, cmap='Reds', annot=True, fmt="1.3f", square=True)

plt.show()
