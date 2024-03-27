import pandas as pd
import numpy as np

def filter_condition(row):
    # 正の値と負の値のカウント
    positive_count = np.sum(row > 0)
    negative_count = np.sum(row < 0)
    
    # 正が1つもない、または負が1つもない条件を満たすかチェック
    return (positive_count == 0) or (negative_count == 0)

df = pd.read_csv('A5SS.csv', sep='\t', skiprows=1, header=None)

# 最初の5列を基に重複を検出
cols_to_check = [0, 1, 2, 3, 4, 5]

# 重複を含む行を検出
duplicates = df.duplicated(subset=cols_to_check, keep=False)

# 重複がある場合の処理
if duplicates.any():
    # 各重複グループに対する処理
    for _, group in df[duplicates].groupby(cols_to_check):
        # グループ内のインデックスを取得
        indexes = group.index.tolist()
        # 元の行に重複の値を追加
        for idx in indexes[1:]:
            # 追加する列を決定（7列目以降に空いている最初の列を探す、なければ新たに作成）
            col_to_add = df.columns.where(df.loc[indexes[0]].isnull()).min()
            if pd.isna(col_to_add):  # 必要な列が存在しない場合は新しい列を追加
                col_to_add = len(df.columns)
            # 値を追加
            df.at[indexes[0], col_to_add] = df.at[idx, 6]
            # 重複した行を削除
            df.drop(index=idx, inplace=True)

multi = df[pd.to_numeric(df[7], errors='coerce').notnull()]

multi.iloc[:, 6:10] = multi.iloc[:, 6:10].apply(pd.to_numeric, errors='coerce')

# 6列目から9列目の行に対してフィルタリング条件を適用
filtered_df = multi.iloc[:, 6:10].apply(filter_condition, axis=1)

# 条件に一致する行のみを含むデータフレームを抽出
result_df = multi[filtered_df]

for index, row in result_df.iterrows():
    if row[5] == "-":
        # 2列目と3列目の値を+1
        result_df.at[index, 1] += 1
        result_df.at[index, 2] += 1
    elif row[5] == "+":
        # 4列目の値を+1
        result_df.at[index, 3] += 1

result_df.to_csv('A5SS_dup3.tsv', sep='\t', index=False, header=None)
