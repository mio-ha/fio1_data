import pandas as pd
import csv

#event = ["A5SS", "A3SS", "RI", "SE"]

event=["SE"]

for ev in event:
	# 最大列数を格納する変数
	max_cols = 0
	input_path = 'combined_{}4.tsv'.format(ev)
	output_path = 'combined_{}4.1.tsv'.format(ev)

	# 最大列数を見つける
	with open(input_path, 'r', encoding='utf-8') as file:
		for line in file:
		    cols = line.count('\t') + 1  # タブで列数をカウント
		    max_cols = max(max_cols, cols)

	# 不足している列を補ってファイルを書き直す
	with open(input_path, 'r', encoding='utf-8') as infile, open(output_path, 'w', newline='', encoding='utf-8') as outfile:
		reader = csv.reader(infile, delimiter='\t')
		writer = csv.writer(outfile, delimiter='\t')
		
		for row in reader:
		    # 不足分の列を追加
		    row += [''] * (max_cols - len(row))
		    writer.writerow(row)


	# TSVファイルの読み込み
	df = pd.read_csv(output_path, sep='\t', header=None)

	# 1-3列目をキーとして重複を除去
	df_unique = df.sort_values(by=[df.columns[0], df.columns[1], df.columns[2]]).drop_duplicates(subset=[df.columns[0], df.columns[1], df.columns[2]])

	# 6列目の値に基づいて負の行と正の行を分ける
	if ev == "RI" or ev == "SE":
		df_negative = df_unique[df_unique[df.columns[5]] < 0]
		df_positive = df_unique[df_unique[df.columns[5]] > 0]
	else:
		df_negative = df_unique[df_unique[df.columns[6]] < 0]
		df_positive = df_unique[df_unique[df.columns[6]] > 0]

	# 結果を異なるTSVファイルに保存
	df_negative.to_csv('sorted_unique_{}4_minus.tsv'.format(ev), sep='\t', index=False, header=False)
	df_positive.to_csv('sorted_unique_{}4_plus.tsv'.format(ev), sep='\t', index=False, header=False)
