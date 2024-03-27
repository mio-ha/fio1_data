import pandas as pd
import statistics as st

all = pd.read_table("edgeR/edgeR_result_splicing_Parker_sig.csv")

fil_all = all[all["FDR"] < 0.05]

sam2 = [[4, 5], [2, 4], [4, 6], [3, 4], [1, 2], [5, 6]]
sam4 = [[4, 6], [1, 4], [1, 5], [2, 5], [2, 4], [1, 3]]

sum_rep4 = 0
sum_rep2 = 0
list_rep4 = []
list_rep2 = []

for i in range(1):
	#rep4 = pd.read_table("edgeR/4rep/edgeR_result_splicing_Parker_4rep_{}_{}.txt".format(sam4[i][0], sam4[i][1]))
	#rep2 = pd.read_table("edgeR/2rep/edgeR_result_splicing_Parker_2rep_{}_{}.txt".format(sam2[i][0], sam2[i][1]))
	#rep4 = pd.read_table("edgeR/edgeR_result_splicing_Parker_667_{}.txt".format(i+1), low_memory=False)
	#rep2 = pd.read_table("edgeR/edgeR_result_splicing_Parker_333_{}.txt".format(i+1), low_memory=False)
	rep4 = pd.read_table("edgeR/edgeR_result_splicing_Parker_100bp.txt", low_memory=False)
	rep2 = pd.read_table("edgeR/edgeR_result_splicing_Parker_50bp.txt", low_memory=False)
	fil_rep4 = rep4[rep4["FDR"] < 0.05]
	fil_rep2 = rep2[rep2["FDR"] < 0.05]
	merge4 = pd.merge(fil_all, fil_rep4, on=["Geneid", "Start"], how='inner')
	merge2 = pd.merge(fil_all, fil_rep2, on=["Geneid", "Start"], how='inner')

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
#print("ave. 4rep", ave_rep4/len(fil_all))
#print("ave. 2rep", ave_rep2/len(fil_all))
print("ave. 667", st.mean(list_rep4)/len(fil_all))
print("ave. 333", st.mean(list_rep2)/len(fil_all))
print("med. 667", st.median(list_rep4)/len(fil_all))
print("med. 333", st.median(list_rep2)/len(fil_all))
'''
