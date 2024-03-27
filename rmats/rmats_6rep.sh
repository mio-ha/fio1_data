for i in $(seq 2 6); do

echo "/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081485_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081486_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081487_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081488_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081489_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081490_333Aligned.sortedByCoord.out.bam," > col0_Parker_333_${i}_samples.txt
echo "/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081509_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081510_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081511_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081512_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081513_333Aligned.sortedByCoord.out.bam,/media/miyokawa/8TB-Data3/test26/333_${i}/ERR9081514_333Aligned.sortedByCoord.out.bam," > fio1_Parker_333_${i}_samples.txt

python /home/miyokawa/ドキュメント/ngs/software/rmats_turbo_v4_1_2/rmats.py \
--b1 col0_Parker_333_${i}_samples.txt --b2 fio1_Parker_333_${i}_samples.txt \
--gtf /home/miyokawa/ドキュメント/ngs/araport/Araport11_GTF_genes_transposons.Apr2023_nochr.gtf \
-t paired --variable-read-length --readLength 150 --nthread 4 --cstat 0.1 \
--od rmats_Parker_333_${i} \
--tmp rmats_tmp_Parker_333_${i}

done
