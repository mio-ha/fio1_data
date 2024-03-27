#!/bin/bash

# 優先度を関連付ける
declare -A files
files=(
  ["DEXSeq/DEXSeq_Parker20c_sorted.bed"]=4
  ["edgeR/edgeR_result_splicing_Parker_sorted.bed"]=3
  ["rmats/ALL.MATS.JCEC_fix_Parker20c.bed"]=1
  ["SUPPA/diffSplice/ara_diffSplice_Parker20c_fix.bed"]=2
)

declare -A aliases
aliases=(
  ["DEXSeq/DEXSeq_Parker20c_sorted.bed"]="DEXSeq"
  ["edgeR/edgeR_result_splicing_Parker_sorted.bed"]="edgeR"
  ["rmats/ALL.MATS.JCEC_fix_Parker20c.bed"]="rMATS"
  ["SUPPA/diffSplice/ara_diffSplice_Parker20c_fix.bed"]="SUPPA"
)

# ファイルの名前を配列に追加
file_names=("${!files[@]}")

# 組み合わせを取得してコマンドを実行
for ((i=0; i<${#file_names[@]}; i++)); do
  for ((j=i+1; j<${#file_names[@]}; j++)); do
    # 優先度に基づいてソート
    if (( ${files[${file_names[i]}]} > ${files[${file_names[j]}]} )); then
      bedtools intersect -wa -a "${file_names[i]}" -b "${file_names[j]}" -sorted > "overlap/overlap_reverse/overlap_fix_Parker_${aliases[${file_names[i]}]}_${aliases[${file_names[j]}]}.bed"
    else
      bedtools intersect -wa -a "${file_names[j]}" -b "${file_names[i]}" -sorted > "overlap/overlap_reverse/overlap_fix_Parker_${aliases[${file_names[j]}]}_${aliases[${file_names[i]}]}.bed"
    fi
  done
done

