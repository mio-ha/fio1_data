import pandas as pd

batch = pd.read_table("DEXSeq_all_norem_batch_result.txt", low_memory=False)
no_batch = pd.read_table("DEXSeq_all_rem_batch_result.txt", low_memory=False)

batch = batch[batch["padj"] < 0.05]
no_batch = no_batch[no_batch["padj"] < 0.05]

batch["name"] = batch["groupID"] + batch["featureID"]
no_batch["name"] = no_batch["groupID"] + no_batch["featureID"]

merge = pd.merge(batch, no_batch, on="name", how="inner")
print(merge)
