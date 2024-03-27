import pandas as pd
import statistics as st

all = pd.read_table("DEXSeq/DEXSeq_Parker_sig_abs.csv")

fil_all = all[all["padj"] < 0.05]

sam2 = [[4, 5], [2, 4], [4, 6], [3, 4], [1, 2], [5, 6]]
sam4 = [[4, 6], [1, 4], [1, 5], [2, 5], [2, 4], [1, 3]]

sum_rep4 = 0
sum_rep2 = 0
list_rep4 = []
list_rep2 = []

for i in range(1):
	#rep4 = pd.read_table("DEXSeq/4rep/DEXSeq_Parker_4rep_{}_{}_result.txt".format(sam4[i][0], sam4[i][1]), low_memory=False)
	#rep2 = pd.read_table("DEXSeq/2rep/DEXSeq_Parker_2rep_{}_{}_result.txt".format(sam2[i][0], sam2[i][1]), low_memory=False)
	#rep4 = pd.read_table("DEXSeq/DEXSeq_Parker_667_{}_result.txt".format(i+1), low_memory=False)
	#rep2 = pd.read_table("DEXSeq/DEXSeq_Parker_333_{}_result.txt".format(i+1), low_memory=False)
	rep4 = pd.read_table("DEXSeq/DEXSeq_Parker_100bp_result.txt", low_memory=False)
	rep2 = pd.read_table("DEXSeq/DEXSeq_Parker_50bp_result.txt", low_memory=False)
	fil_rep4 = rep4[rep4["padj"] < 0.05]
	fil_rep2 = rep2[rep2["padj"] < 0.05]
	merge4 = pd.merge(fil_all, fil_rep4, on=["groupID", "featureID"], how='inner')
	merge2 = pd.merge(fil_all, fil_rep2, on=["groupID", "featureID"], how='inner')

	#sum_rep4 += len(merge4)
	#sum_rep2 += len(merge2)
	list_rep4.append(len(merge4))
	list_rep2.append(len(merge2))
	#print("4rep", len(merge4)/len(fil_all))
	#print("2rep", len(merge2)/len(fil_all))
	#print("667", len(merge4)/len(fil_all))
	#print("333", len(merge2)/len(fil_all))
	print("100bp", len(merge4)/len(fil_all))
	print("50bp", len(merge2)/len(fil_all))
	
'''
print("ave. 667", st.mean(list_rep4)/len(fil_all))
print("ave. 333", st.mean(list_rep2)/len(fil_all))
print("med. 667", st.median(list_rep4)/len(fil_all))
print("med. 333", st.median(list_rep2)/len(fil_all))
'''
