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
# 差分の結果を表示
echo ${diff[@]}
echo "1,$((diff[0]+1)),$((diff[1]+1)),$((diff[2]+1)),$((diff[3]+1))"

cut -f 1,$((diff[0]+1)),$((diff[1]+1)),$((diff[2]+1)),$((diff[3]+1)) psi/Parker20c_event_psi_col0.txt > 4rep/control_${i}_${j}.psi
cut -f 1,$((diff[0]+1)),$((diff[1]+1)),$((diff[2]+1)),$((diff[3]+1)) tpm/Parker20c_tpm_col0.txt > 4rep/control_${i}_${j}.tpm

cut -f 1,$((diff[0]+1)),$((diff[1]+1)),$((diff[2]+1)),$((diff[3]+1)) psi/Parker20c_event_psi_fio1.txt > 4rep/treat_${i}_${j}.psi
cut -f 1,$((diff[0]+1)),$((diff[1]+1)),$((diff[2]+1)),$((diff[3]+1)) tpm/Parker20c_tpm_fio1.txt > 4rep/treat_${i}_${j}.tpm

python /home/miyokawa/ドキュメント/ngs/software/SUPPA-2.3/suppa.py diffSplice \
-m empirical -gc \
-i ioe/ara.all.events.ioe \
--save_tpm_events \
-p 4rep/treat_${i}_${j}.psi 4rep/control_${i}_${j}.psi \
-e 4rep/treat_${i}_${j}.tpm 4rep/control_${i}_${j}.tpm \
-o 4rep/ara_diffSplice_Parker20c_4rep_${i}_${j}

done
done
