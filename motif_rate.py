import pandas as pd

df = pd.read_table("SUPPA/ioe/ara.all.events2.ioe", header=None, skiprows=1)

#A3
target = df[df.iloc[:, 3] == "A3"]
strand = target[target.iloc[:, 9] == "+"]
A3_5 = strand.iloc[:, 5]
A3_3_A1 = strand.iloc[:, 6]
A3_3_A2 = strand.iloc[:, 8]

strand = target[target.iloc[:, 9] == "-"]
A3_5 = pd.concat([A3_5, strand.iloc[:, 6]], axis=0, ignore_index=True).astype(int)
A3_3_A1 = pd.concat([A3_3_A1, strand.iloc[:, 5]], axis=0, ignore_index=True).astype(int)
A3_3_A2 = pd.concat([A3_3_A2, strand.iloc[:, 7]], axis=0, ignore_index=True).astype(int)

#A5
target = df[df.iloc[:, 3] == "A5"]
strand = target[target.iloc[:, 9] == "+"]
A5_5_A1 = strand.iloc[:, 7]
A5_5_A2 = strand.iloc[:, 5]
A5_3 = strand.iloc[:, 6]

strand = target[target.iloc[:, 9] == "-"]
A5_5_A1 = pd.concat([A5_5_A1, strand.iloc[:, 6]], axis=0, ignore_index=True).astype(int)
A5_5_A2 = pd.concat([A5_5_A2, strand.iloc[:, 8]], axis=0, ignore_index=True).astype(int)
A5_3 = pd.concat([A5_3, strand.iloc[:, 5]], axis=0, ignore_index=True).astype(int)

#RI
target = df[df.iloc[:, 3] == "RI"]
strand = target[target.iloc[:, 9] == "+"]
RI_5 = strand.iloc[:, 5]
RI_3 = strand.iloc[:, 6]

strand = target[target.iloc[:, 9] == "-"]
RI_5 = pd.concat([RI_5, strand.iloc[:, 6]], axis=0, ignore_index=True).astype(int)
RI_3 = pd.concat([RI_3, strand.iloc[:, 5]], axis=0, ignore_index=True).astype(int)

#SE
target = df[df.iloc[:, 3] == "SE"]
strand = target[target.iloc[:, 9] == "+"]
SE_5 = strand.iloc[:, 5]
SE_3 = strand.iloc[:, 6]

strand = target[target.iloc[:, 9] == "-"]
SE_5 = pd.concat([RI_5, strand.iloc[:, 6]], axis=0, ignore_index=True).astype(int)
SE_3 = pd.concat([RI_3, strand.iloc[:, 5]], axis=0, ignore_index=True).astype(int)

print(RI_3, RI_5)
