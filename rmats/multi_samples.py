# ベースとなるテキストフォーマット（プレースホルダーを使用）
base_text_format = "/media/miyokawa/8TB-Data3/STAR_sample/ERR{}_100bpAligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/STAR_sample/ERR{}_100bpAligned.sortedByCoord.out.bam"

# 組み合わせリスト
pairs = [(3, 5), (2, 6), (3, 6), (2, 5), (1, 6), (4, 5)]

# ファイル名のフォーマット
file_name_format = "Col0_Parker_100bp_2rep_samples_{}.txt"

pairl=[0,0]

for i, pair in enumerate(pairs, start=1):
    pairl[0] = 9081484 + pair[0] #9081484
    pairl[1] = 9081484 + pair[1] #9081508

    # 各組み合わせに基づいてテキストを生成（.formatを使用）
    text = base_text_format.format(pairl[0], pairl[1])
    
    # ファイル名を生成
    file_name = file_name_format.format(i)
    
    # テキスト内容をファイルに書き込む
    with open(file_name, 'w') as file:
        file.write(text + '\n')

