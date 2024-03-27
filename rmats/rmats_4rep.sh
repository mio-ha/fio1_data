for i in $(seq 1 6); do
for j in $(seq 1 6); do

if [[ $i -ge $j ]]; then
	continue
fi

A=($(seq 1 6))
# 配列Bを定義
B=($i $j)
# 配列Aと配列Bの差分を求める
diff=()
for element in "${A[@]}"; do
  found=0
  for compare in "${B[@]}"; do
    if [[ "$element" == "$compare" ]]; then
      found=1
      break
    fi
  done
  if [[ $found -eq 0 ]]; then
    diff+=("$element")
  fi
done

echo "../tech_bam/Parker_Col-0_20c_${diff[0]}.bam,../tech_bam/Parker_Col-0_20c_${diff[1]}.bam,../tech_bam/Parker_Col-0_20c_${diff[2]}.bam,../tech_bam/Parker_Col-0_20c_${diff[3]}.bam" > col0_Parker_4rep_${i}_${j}_samples.txt
echo "../tech_bam/Parker_fio1_20c_${diff[0]}.bam,../tech_bam/Parker_fio1_20c_${diff[1]}.bam,../tech_bam/Parker_fio1_20c_${diff[2]}.bam,../tech_bam/Parker_fio1_20c_${diff[3]}.bam" > fio1_Parker_4rep_${i}_${j}_samples.txt

python /home/miyokawa/ドキュメント/ngs/software/rmats_turbo_v4_1_2/rmats.py \
--b1 col0_Parker_4rep_${i}_${j}_samples.txt --b2 fio1_Parker_4rep_${i}_${j}_samples.txt \
--gtf /home/miyokawa/ドキュメント/ngs/araport/Araport11_GTF_genes_transposons.Apr2023_nochr.gtf \
-t paired --variable-read-length --readLength 150 --nthread 4 --cstat 0.1 \
--od rmats_Parker_4rep_${i}_${j} \
--tmp rmats_tmp_Parker_4rep_${i}_${j}

done
done
