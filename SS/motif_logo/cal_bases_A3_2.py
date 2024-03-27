from Bio import SeqIO
import pandas as pd
import logomaker

def extract_sequences(fasta_file, bed_file):
    sequences = {'start': [], 'end1': [], 'end2': []}
    for record in SeqIO.parse(fasta_file, "fasta"):
        with open(bed_file, 'r') as bed:
            for line in bed:
                parts = line.strip().split()
                chromosome, strand = parts[0], parts[5]
                start, end1, end2 = int(parts[1]), int(parts[2]), int(parts[3])
                
                if record.id != chromosome:
                    continue

                # 配列の抽出
                extract_seq = lambda pos: record.seq[max(0, pos-4):pos+5].upper()
                extract_seq_rev = lambda pos: record.seq[max(0, pos-4):pos+5].reverse_complement().upper()
                
                # 鎖に応じて配列を抽出
                sequences['start'].append(extract_seq_rev(start-2) if strand == '-' else extract_seq(start))
                sequences['end1'].append(extract_seq_rev(end1) if strand == '-' else extract_seq(end1-2))
                sequences['end2'].append(extract_seq_rev(end2) if strand == '-' else extract_seq(end2-2))
                
    return sequences

def create_info_matrix(sequences):
    info_matrices = {}
    for key, seq_list in sequences.items():
        nt_counts = pd.DataFrame(0, index=['A', 'C', 'G', 'T'], columns=range(-4, 5))
        for seq in seq_list:
            for i, nt in enumerate(seq):
                if nt in nt_counts.index:
                    nt_counts.loc[nt, i-4] += 1
        nt_totals = nt_counts.sum()
        info_matrices[key] = nt_counts.div(nt_totals, axis=1)
        #info_matrices[key] = nt_counts
    return info_matrices

def generate_sequence_logo(info_matrices, outfile):
    for key, info_matrix in info_matrices.items():
        if key == 'start':
            position = ["−4", "−3", "−2", "−1", "+1", "+2", "+3", "+4", "+5"]
        else:
            position = ["−5", "−4", "−3", "−2", "−1", "+1", "+2", "+3", "+4"]
        logo = logomaker.Logo(info_matrix.T, fade_below=0.5)
        logo.ax.tick_params(axis='both', labelsize=16)  # X軸とY軸の両方の軸目盛り
        logo.ax.set_xticks(range(-4, len(position)-4))  # X軸の目盛り位置を設定
        logo.ax.set_xticklabels(position) 
        logo.ax.set_ylabel("Probability", fontsize=16)
        #logo.style_xticks(rotation=90)
        #logo.ax.set_title(f"{key} Position")
        logo.fig.savefig(f"{key}{outfile}")

fasta_file = "/home/miyokawa/ドキュメント/ngs/araport/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"  # FASTAファイルのパス
bed_file = ["SS/sorted_unique_A3SS4_plus.tsv", "SS/sorted_unique_A3SS4_minus.tsv"]  # BEDファイルのパス
out_file = ["_position_A3_plus_logo4.png", "_position_A3_minus_logo4.png"]

for i in range(2):
	# 開始位置1、開始位置2、終了位置の周辺領域から配列を抽出
	sequences = extract_sequences(fasta_file, bed_file[i])

	# 各領域の配列から情報行列を作成
	info_matrices = create_info_matrix(sequences)
	
	#TをUに変換
	for key, frame in info_matrices.items():
	    frame.rename(index={"T": "U"}, inplace=True)
	    
	# シークエンスロゴを生成してファイルに保存
	generate_sequence_logo(info_matrices, out_file[i])

