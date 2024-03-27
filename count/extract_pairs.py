import pandas as pd

# ファイルの読み込みパス
file_path = 'dex_counts_down_flat2.txt'

# 抽出する列のインデックスを計算する関数
def calculate_column_indices_corrected(pair):
    # 7列目が1、12列目が6なので、それに合わせて調整
    # ここでは、7から始まる範囲のみを考慮する
    return [5 + pair[0], 5 + pair[1], 11 + pair[0], 11 + pair[1]]

# ファイルを読み込む
data = pd.read_csv(file_path, sep='\t', skiprows=1)

# 各組み合わせに対して列を追加
pairs = [(3, 5), (2, 6), (3, 6), (2, 5), (1, 6), (4, 5)]
i = 0
for pair in pairs:
	i = i + 1
	# 初期の列（1-6列目）
	selected_columns_corrected = list(range(6))
	selected_columns_corrected.extend(calculate_column_indices_corrected(pair))
	print(selected_columns_corrected)
	# 対応する列を抽出
	extracted_data_corrected = data.iloc[:, selected_columns_corrected]

	# 保存
	extracted_data_corrected.to_csv("dex_counts_down_flat_{}.txt".format(i), sep="\t", index=False)

