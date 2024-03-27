for i in $(seq 7 12)
do
for j in $(seq 7 12)
do
if [ $i -ge $j ]; then
continue
fi

cut -f 1-6,${i},${j},$((i+6)),$((j+6)) dex_counts_Parker20c.txt > dex_counts_Parker20c_2rep_${i}_${j}.txt

done
done
