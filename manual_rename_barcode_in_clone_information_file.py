import os
import re
import pandas as pd
import sys

input_barcode = sys.argv[1]
input_clone = sys.argv[2]
out_file = sys.argv[3]

# obs_file="new_obs_file.csv"
# df_obs=pd.read_csv(obs_file,index_col=0)
df_barcodes=pd.read_table(input_barcode)
df_barcodes["rename_obs_barcodes"]=[ "{}_{}".format(a.split("-")[0],a.split("_")[1]) for a in df_barcodes["x"].values]
df_barcodes["rename_clone_barcodes"]=[ a.split("_")[1] for a in df_barcodes["rename_obs_barcodes"].values]
df_clone=pd.read_csv(input_clone)


# obs_merge=df_obs.merge(df_barcodes,left_on="barcodes",right_on="rename_obs_barcodes",how="inner")
# obs_merge.drop(columns=["barcodes","rename_obs_barcodes","rename_clone_barcodes"],inplace=True)
# obs_merge["barcodes"]=obs_merge["x"]
# obs_merge.set_index(keys=["x"],drop=True,inplace=True)
# obs_merge.to_csv("obs_merge.csv",index=True)

clone_merge=df_barcodes.merge(df_clone,left_on="rename_clone_barcodes",right_on="sample_id",how="inner")
clone_merge.drop(columns=["rename_obs_barcodes","rename_clone_barcodes","sample_id"],inplace=True)
clone_merge.columns=["sample_id","clone_id"]
print(clone_merge)
clone_merge.to_csv(out_file,index=False)

