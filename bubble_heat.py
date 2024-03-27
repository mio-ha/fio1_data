import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def text_color_based_on_background(color, threshold=0.5):
    if color < threshold:
        return 'black'
    else:
        return 'white'

#pipeline = ["DEXSeq", "edgeR", "rMATS", "SUPPA"]
pipeline = ["Cai", "Parker", "Sun", "Wang"]

# データのロード
for pipe in pipeline:
    heatmap_data = pd.read_csv('overlap_table_method_{}_J.csv'.format(pipe), sep='\t', index_col=0)
    bubble_data = pd.read_csv('overlap_table_method_count_{}3.csv'.format(pipe), sep='\t', index_col=0)

    heatmap_data = heatmap_data.iloc[::-1]
    bubble_data = bubble_data.iloc[::-1]

    # NumPy配列に変換
    heatmap_values = heatmap_data.values
    bubble_values = bubble_data.values

    # ヒートマップの値を0-1の範囲にスケーリング（カラーマップを適用するため）
    heatmap_values_scaled = heatmap_values#(heatmap_values - np.min(heatmap_values)) / (np.max(heatmap_values) - np.min(heatmap_values))

    # カラーマップを定義
    cmap1 = plt.get_cmap('Blues')
    cmap2 = plt.get_cmap('Reds')

    # 最大バブル値を定義
    max_bubble_value = 1654 #40.67 #1654 #820
    bubble_scale = 75**2 #75**3
    bubble_threshold = 40

    # バブルサイズのスケーリング
    bubble_sizes = np.sqrt(bubble_scale * (bubble_values / max_bubble_value))
    #bubble_sizes = np.sqrt(np.sqrt(bubble_scale * (bubble_values / max_bubble_value)))

    # グラフを描画（正方形にするためにfigsizeを調整）
    fig, ax = plt.subplots(figsize=[6.5, 6])
    plt.rcParams["font.family"] = "Liberation Sans"
    plt.rcParams["font.size"] = 18

    # 目盛り線を消去
    ax.grid(False)

    # バブルチャートを描画
    for i in range(bubble_values.shape[0]):
        for j in range(bubble_values.shape[1]):
        # バブルのサイズが0より大きい場合のみ描画
            if bubble_values[i, j] > 0:
                # ヒートマップの値に基づいてバブルの色を設定
                if i+j < bubble_values.shape[0]:  # 左下部分のとき
                    bubble_color = cmap1(heatmap_values_scaled[i, j])
                else:  # 右上部分のとき
                    bubble_color = cmap2(heatmap_values_scaled[i, j])
                # plot関数を使用してバブルを円形に描画
                ax.plot(j, i, marker='o', markersize=bubble_sizes[i, j], linestyle='', markeredgecolor='white', markerfacecolor=bubble_color)
                # ヒートマップの値をテキストで表示
                if bubble_sizes[i, j] < bubble_threshold:
                    ax.text(j, i+0.3, f"{heatmap_values[i, j]:.2f}", ha='center', va='center', color='black')
                else:
                     ax.text(j, i, f"{heatmap_values[i, j]:.2f}", ha='center', va='center', color=text_color_based_on_background(heatmap_values_scaled[i, j]))

    # ラベルを追加
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top') 
    ax.set_xticks(np.arange(len(heatmap_data.columns)))
    ax.set_yticks(np.arange(len(heatmap_data.index)))
    
    heatmap_data.columns = heatmap_data.columns.str.replace("\\n", "\n")
    heatmap_data.index = heatmap_data.index.str.replace("\\n", "\n")
    ax.set_xticklabels(heatmap_data.columns, fontsize=14)
    ax.set_yticklabels(heatmap_data.index, fontsize=14)

    # 項目と項目の中間に線を入れる（線の色を灰色にする）
    ax.set_xticks(np.arange(-.5, len(heatmap_data.columns), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(heatmap_data.index), 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    plt.subplots_adjust(left=0.155) #0.155
    plt.savefig("{}_method_bubble_J.pdf".format(pipe), dpi=300)
    #plt.show()
