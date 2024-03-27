for i in $(seq 1 6); do
for j in $(seq 1 6); do

if [[ $i -ge $j ]]; then
	continue
fi

cut -f 1,$((i+1)),$((j+1)) psi/Parker20c_event_psi_col0.txt > 2rep/control_${i}_${j}.psi
cut -f 1,$((i+1)),$((j+1)) tpm/Parker20c_tpm_col0.txt > 2rep/control_${i}_${j}.tpm

cut -f 1,$((i+1)),$((j+1)) psi/Parker20c_event_psi_fio1.txt > 2rep/treat_${i}_${j}.psi
cut -f 1,$((i+1)),$((j+1)) tpm/Parker20c_tpm_fio1.txt > 2rep/treat_${i}_${j}.tpm

python /home/miyokawa/ドキュメント/ngs/software/SUPPA-2.3/suppa.py diffSplice \
-m empirical -gc \
-i ioe/ara.all.events.ioe \
--save_tpm_events \
-p 2rep/treat_${i}_${j}.psi 2rep/control_${i}_${j}.psi \
-e 2rep/treat_${i}_${j}.tpm 2rep/control_${i}_${j}.tpm \
-o 2rep/ara_diffSplice_Parker20c_2rep_${i}_${j}

done
done
