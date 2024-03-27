import pandas as pd
import io

# ファイルを行ごとに読み込みます
with open("rmats_Wang/ALL.MATS.JCEC.txt", 'r') as f:
    lines = f.readlines()

# データフレームを格納するための空のリストを作成します
dfs = []
current_table = []

# 各行をチェックします
for line in lines:
    # 新しいテーブルが始まるかどうかをチェックします
    if line.startswith("ID"):
        # 現在のテーブルが空でない場合、それをデータフレームに変換します
        if current_table:
            df = pd.read_csv(io.StringIO("\n".join(current_table)), sep="\t")
            dfs.append(df)
        # 新しいテーブルを開始します
        current_table = [line.strip()]
    else:
        current_table.append(line.strip())

# 最後のテーブルをデータフレームに変換します
if current_table:
    df = pd.read_csv(io.StringIO("\n".join(current_table)), sep="\t")
    dfs.append(df)
'''
# 最初の5つのデータフレームを表示します
for i, df in enumerate(dfs[:5]):
    print(f"DataFrame {i+1}:")
    print(df.head())
    print("\n---\n")
'''
for i in range(len(dfs[0])):
    if dfs[0].loc[i, "strand"] == "+":
        dfs[0]["Start"] = dfs[0]["longExonStart_0base"]
        dfs[0]["End"] = dfs[0]["shortES"]
    else:
        dfs[0]["Start"] = dfs[0]["shortEE"]
        dfs[0]["End"] = dfs[0]["longExonEnd"]
        
for i in range(len(dfs[1])):
    if dfs[1].loc[i, "strand"] == "+":
        dfs[1]["Start"] = dfs[1]["shortEE"]
        dfs[1]["End"] = dfs[1]["longExonEnd"]
    else:
        dfs[1]["Start"] = dfs[1]["longExonStart_0base"]
        dfs[1]["End"] = dfs[1]["shortES"]

dfs[2]["Start"] = dfs[2]["1stExonStart_0base"]
dfs[2]["End"] = dfs[2]["2ndExonEnd"]

dfs[3]["Start"] = dfs[3]["upstreamEE"]
dfs[3]["End"] = dfs[3]["downstreamES"]

dfs[4]["Start"] = dfs[4]["exonStart_0base"]
dfs[4]["End"] = dfs[4]["exonEnd"]

exdfs = dfs

for i in range(5):
    dfs[i]["abs"] = abs(dfs[i]["IncLevelDifference"])
    dfs[i]["chr"] = dfs[i]["chr"].str.replace("chr", "")
    exdfs[i] = dfs[i][["chr", "Start", "End", "GeneID", "FDR", "strand", "IncLevelDifference", "abs"]]

out = pd.concat(exdfs[0:6])
out.sort_values(["chr", "Start"], inplace=True)
out.to_csv("ALL.MATS.JCEC_fix_Wang.bed", sep="\t", index=False, header=False)
