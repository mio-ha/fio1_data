import pandas as pd
import sys

# ファイルをデータフレームとして読み込みます
dff = pd.read_csv("diffSplice/ara_diffSplice_Wang_st.csv", sep="\t", header=None, skiprows=1)

# 'event'列の値に基づいてデータフレームをグループ分けします
grouped = dff.groupby(1)

# グループごとにデータフレームを抽出し、それを辞書に格納します
dfs_by_event = {event: group.reset_index(drop=True) for event, group in grouped}

# 各データフレームについて、ハイフンでつながっている数字を別の列に分割します
for event in dfs_by_event:
    # 各列についてチェックします
    for col in dfs_by_event[event].columns:
        # データの型が文字列の列のみを対象にします
        if dfs_by_event[event][col].dtype == object:
            # 列のデータがハイフンでつながっている数字であるかどうかをチェックします
            # さらに、ハイフンでつながっている数字が完全に列のデータを構成しているかどうかもチェックします
            if dfs_by_event[event][col].str.match('^[0-9]+-[0-9]+$').any():
                # ハイフンでつながっている数字を別の列に分割します
                dfs_by_event[event][[f'{col}_1', f'{col}_2']] = dfs_by_event[event][col].str.split('-', expand=True)
                # 元の列を削除します
                dfs_by_event[event].drop(col, axis=1, inplace=True)
'''
# 再度、識別記号ごとの最初のデータフレームを表示します
for event, df in dfs_by_event.items():
    print(f"DataFrame for event '{event}':")
    print(df.head())
    print("\n---\n")
'''
#pd.set_option("display.max_columns", None)

dfs_by_event["A3"]["p-value"] = dfs_by_event["A3"][7]
dfs_by_event["A5"]["p-value"] = dfs_by_event["A5"][7]
dfs_by_event["AF"]["p-value"] = dfs_by_event["AF"][9]
dfs_by_event["AL"]["p-value"] = dfs_by_event["AL"][9]
dfs_by_event["MX"]["p-value"] = dfs_by_event["MX"][9]
dfs_by_event["RI"]["p-value"] = dfs_by_event["RI"][8]
dfs_by_event["SE"]["p-value"] = dfs_by_event["SE"][7]
dfs_by_event["A3"]["strand"] = dfs_by_event["A3"][5]
dfs_by_event["A5"]["strand"] = dfs_by_event["A5"][5]
dfs_by_event["AF"]["strand"] = dfs_by_event["AF"][7]
dfs_by_event["AL"]["strand"] = dfs_by_event["AL"][7]
dfs_by_event["MX"]["strand"] = dfs_by_event["MX"][7]
dfs_by_event["RI"]["strand"] = dfs_by_event["RI"][6]
dfs_by_event["SE"]["strand"] = dfs_by_event["SE"][5]
dfs_by_event["A3"]["psi"] = dfs_by_event["A3"][6]
dfs_by_event["A5"]["psi"] = dfs_by_event["A5"][6]
dfs_by_event["AF"]["psi"] = dfs_by_event["AF"][8]
dfs_by_event["AL"]["psi"] = dfs_by_event["AL"][8]
dfs_by_event["MX"]["psi"] = dfs_by_event["MX"][8]
dfs_by_event["RI"]["psi"] = dfs_by_event["RI"][7]
dfs_by_event["SE"]["psi"] = dfs_by_event["SE"][6]

for event in dfs_by_event:
    dfs_by_event[event] = dfs_by_event[event][dfs_by_event[event]["p-value"].astype(float) < 0.05].reset_index(drop=True)
    dfs_by_event[event]["Start"] = None
    dfs_by_event[event]["End"] = None
    dfs_by_event[event]["abs"] = abs(dfs_by_event[event]["psi"].astype(float))
    #print(len(dfs_by_event[event]))

#sys.exit()

for i in range(len(dfs_by_event["A3"])):
    if dfs_by_event["A3"].loc[i, "strand"] == "+":
        dfs_by_event["A3"].loc[i, "Start"] = dfs_by_event["A3"].loc[i, "3_2"]
        dfs_by_event["A3"].loc[i, "End"] = dfs_by_event["A3"].loc[i, "4_2"]
    else:
        dfs_by_event["A3"].loc[i, "Start"] = dfs_by_event["A3"].loc[i, "4_1"]
        dfs_by_event["A3"].loc[i, "End"] = dfs_by_event["A3"].loc[i, "3_1"]

for i in range(len(dfs_by_event["A5"])):
    if dfs_by_event["A5"].loc[i, "strand"] == "+":
        dfs_by_event["A5"].loc[i, "Start"] = dfs_by_event["A5"].loc[i, "4_1"]
        dfs_by_event["A5"].loc[i, "End"] = dfs_by_event["A5"].loc[i, "3_1"]
    else:
        dfs_by_event["A5"].loc[i, "Start"] = dfs_by_event["A5"].loc[i, "3_2"]
        dfs_by_event["A5"].loc[i, "End"] = dfs_by_event["A5"].loc[i, "4_2"]

for i in range(len(dfs_by_event["AF"])):
    if dfs_by_event["AF"].loc[i, "strand"] == "+":
        dfs_by_event["AF"].loc[i, "Start"] = dfs_by_event["AF"].loc[i, "3_1"]
        dfs_by_event["AF"].loc[i, "End"] = dfs_by_event["AF"].loc[i, "6_1"]
    else:
        dfs_by_event["AF"].loc[i, "Start"] = dfs_by_event["AF"].loc[i, "3_2"]
        dfs_by_event["AF"].loc[i, "End"] = dfs_by_event["AF"].loc[i, "6_1"]
        
for i in range(len(dfs_by_event["AL"])):
    if dfs_by_event["AL"].loc[i, "strand"] == "+":
        dfs_by_event["AL"].loc[i, "Start"] = dfs_by_event["AL"].loc[i, "3_2"]
        dfs_by_event["AL"].loc[i, "End"] = dfs_by_event["AL"].loc[i, "6_1"]
    else:
        dfs_by_event["AL"].loc[i, "Start"] = dfs_by_event["AL"].loc[i, "3_1"]
        dfs_by_event["AL"].loc[i, "End"] = dfs_by_event["AL"].loc[i, "6_1"]
        
for i in range(len(dfs_by_event["MX"])):
    dfs_by_event["MX"].loc[i, "Start"] = dfs_by_event["MX"].loc[i, "3_2"]
    dfs_by_event["MX"].loc[i, "End"] = dfs_by_event["MX"].loc[i, "6_1"]
    
for i in range(len(dfs_by_event["RI"])):
    dfs_by_event["RI"].loc[i, "Start"] = dfs_by_event["RI"].loc[i, "4_1"]
    dfs_by_event["RI"].loc[i, "End"] = dfs_by_event["RI"].loc[i, "4_2"]
    
for i in range(len(dfs_by_event["SE"])):
    dfs_by_event["SE"].loc[i, "Start"] = dfs_by_event["SE"].loc[i, "3_2"]
    dfs_by_event["SE"].loc[i, "End"] = dfs_by_event["SE"].loc[i, "4_1"]

for event in dfs_by_event:
    dfs_by_event[event] = dfs_by_event[event].reindex([2, "Start", "End", 0, "p-value", "strand", "psi", "abs", 1], axis=1)

out = pd.concat(dfs_by_event.values(), keys=dfs_by_event.keys())
out["Start"] = out["Start"].astype(int)
out.sort_values([2, "Start"], inplace=True)
out.to_csv("diffSplice/ara_diffSplice_Wang_fix.bed", sep="\t", index=False, header=False)

