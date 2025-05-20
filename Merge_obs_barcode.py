import pandas as pd
import os
import sys

file_obs=sys.argv[1]
file_barcodes=sys.argv[2]
patient=sys.argv[3]
sample=sys.argv[4]
out_file=sys.argv[5]

df_obs = pd.read_csv(file_obs,header=0,index_col=0)
df_barcode = pd.read_csv(file_barcodes,header=None)

df_barcode.loc[:,"barcodes"]=["{}_{}".format(patient,i) for i in df_barcode.iloc[:,0].tolist()]

df_barcode.columns=["raw","barcodes"]

df_barcode.merge(df_obs,left_on="barcodes",right_on="barcodes",how="inner").loc[1:,"raw"].to_csv(out_file,header=False,index=False)