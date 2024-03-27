for i in $(seq 1 6); do
for j in $(seq 1 6); do

if [[ $i -ge $j ]]; then
	continue
fi

echo "../tech_bam/Parker_Col-0_20c_${i}.bam,../tech_bam/Parker_Col-0_20c_${j}.bam" > col0_Parker_2rep_${i}_${j}_samples.txt
echo "../tech_bam/Parker_fio1_20c_${i}.bam,../tech_bam/Parker_fio1_20c_${j}.bam" > fio1_Parker_2rep_${i}_${j}_samples.txt

python /home/miyokawa/ドキュメント/ngs/software/rmats_turbo_v4_1_2/rmats.py \
--b1 col0_Parker_2rep_${i}_${j}_samples.txt --b2 fio1_Parker_2rep_${i}_${j}_samples.txt \
--gtf /home/miyokawa/ドキュメント/ngs/araport/Araport11_GTF_genes_transposons.Apr2023_nochr.gtf \
-t paired --variable-read-length --readLength 150 --nthread 4 --cstat 0.1 \
--od rmats_Parker_2rep_${i}_${j} \
--tmp rmats_tmp_Parker_2rep_${i}_${j}

done
done
