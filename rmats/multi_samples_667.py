sample = ["Col-0", "fio1"]

for s in sample:
	# ベースとなるテキストフォーマット（プレースホルダーを使用）
	base_text_format = "/media/miyokawa/8TB-Data3/tech_data/tech_bam/test/Parker_{}_20c_{}_667.bam,/media/miyokawa/8TB-Data3/tech_data/tech_bam/test/Parker_{}_20c_{}_667.bam"

	# 組み合わせリスト
	pairs = [(3, 5), (2, 6), (3, 6), (2, 5), (1, 6), (4, 5)]

	# ファイル名のフォーマット
	file_name_format = "{}_Parker_667_2rep_samples_{}.txt"

	for i, pair in enumerate(pairs, start=1):
		# 各組み合わせに基づいてテキストを生成（.formatを使用）
		text = base_text_format.format(s, pair[0], s, pair[1])
		
		# ファイル名を生成
		file_name = file_name_format.format(s, i)
		
		# テキスト内容をファイルに書き込む
		with open(file_name, 'w') as file:
		    file.write(text + '\n')

