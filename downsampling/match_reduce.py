import pandas as pd
import statistics as st

all = pd.read_table("SUPPA/ara_diffSplice_Parker20c.dpsi", header=None, skiprows=1)

fil_all = all[all[2] < 0.05]

sam2 = [[4, 5], [2, 4], [4, 6], [3, 4], [1, 2], [5, 6]]
sam4 = [[4, 6], [1, 4], [1, 5], [2, 5], [2, 4], [1, 3]]

sum_rep4 = 0
sum_rep2 = 0
list_rep4 = []
list_rep2 = []

#for i in range(6):
for i in range(1):
	#rep4 = pd.read_table("SUPPA/4rep/ara_diffSplice_Parker20c_4rep_{}_{}.dpsi.temp.0".format(sam4[i][0], sam4[i][1]), header=None, skiprows=1)
	#rep2 = pd.read_table("SUPPA/2rep/ara_diffSplice_Parker20c_2rep_{}_{}.dpsi.temp.0".format(sam2[i][0], sam2[i][1]), header=None, skiprows=1)
	#rep4 = pd.read_table("SUPPA/667_333/ara_diffSplice_Parker20c_667_{}.dpsi".format(i+1), header=None, skiprows=1)
	#rep2 = pd.read_table("SUPPA/667_333/ara_diffSplice_Parker20c_333_{}.dpsi".format(i+1), header=None, skiprows=1)
	rep4 = pd.read_table("SUPPA/ara_diffSplice_Parker20c_100bp.dpsi", header=None, skiprows=1)
	rep2 = pd.read_table("SUPPA/ara_diffSplice_Parker20c_50bp.dpsi", header=None, skiprows=1)
	fil_rep4 = rep4[rep4[2] < 0.05]
	fil_rep2 = rep2[rep2[2] < 0.05]
	merge4 = pd.merge(fil_all, fil_rep4, on=0, how='inner')
	merge2 = pd.merge(fil_all, fil_rep2, on=0, how='inner')

	list_rep4.append(len(merge4))
	list_rep2.append(len(merge2))
	#print("4rep", len(merge4)/len(fil_all))
	#print("2rep", len(merge2)/len(fil_all))
	#print("667", len(merge4)/len(fil_all))
	#print("333", len(merge2)/len(fil_all))
	print("100bp", len(merge4)/len(fil_all))
	print("50bp", len(merge2)/len(fil_all))
'''	
print("ave. 4rep", st.mean(list_rep4)/len(fil_all))
print("ave. 2rep", st.mean(list_rep2)/len(fil_all))
print("med. 4rep", st.median(list_rep4)/len(fil_all))
print("med. 2rep", st.median(list_rep2)/len(fil_all))
#print("ave. 667", st.mean(list_rep4)/len(fil_all))
#print("ave. 333", st.mean(list_rep2)/len(fil_all))
#print("med. 667", st.median(list_rep4)/len(fil_all))
#print("med. 333", st.median(list_rep2)/len(fil_all))
'''
