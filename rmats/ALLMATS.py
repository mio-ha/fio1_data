import pandas as pd
import os

# 処理するファイルの名前
files = ["A3SS.MATS.JCEC.txt", "A5SS.MATS.JCEC.txt", "MXE.MATS.JCEC.txt", "RI.MATS.JCEC.txt", "SE.MATS.JCEC.txt"]

# ディレクトリ名のパターン
dir_pattern = "rmats_667_100bp_{}"

for i in range(1, 7):
    input_dir = dir_pattern.format(i)
    output_file = os.path.join(input_dir, f"ALL.MATS.JCEC.txt")

    first_file = True
    for file in files:
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path, sep="\t")
        filtered_df = df[df["FDR"] < 0.05]

        # ヘッダー付きでファイルに追記
        with open(output_file, 'a' if not first_file else 'w') as outfile:
            filtered_df.to_csv(outfile, index=False, sep="\t", header=True)
            first_file = False

