for i in $(seq 7 12)
do
for j in $(seq 7 12)
do
if [ $i -ge $j ]; then
continue
fi

cut -f ${i},${j},$((i+6)),$((j+6)) --complement dex_counts_Parker20c.txt > 4rep_matrix/dex_counts_Parker20c_4rep_${i}_${j}.txt

done
done
