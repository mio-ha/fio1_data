import subprocess

study = ["DEXSeq2", "edgeR2", "rMATS2", "SUPPA2.1"]

for i in range(len(study)):
	for j in range(len(study)):
		if j <= i:
			continue
		
		cmd = "bedtools intersect -a ag_logFC_{}_sorted.bed -b ag_logFC_{}_sorted.bed -sorted > ./overlap/overlap_{}_{}_agscore.bed".format(study[i], study[j], study[i], study[j])
		subprocess.call(cmd, shell=True)
		
