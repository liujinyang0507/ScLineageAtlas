import pandas as pd
import numpy as np
import sys

af=sys.argv[1]
gt=sys.argv[2]
cutoff=sys.argv[3]
df_af=pd.read_csv(af,header=0,index_col=0)
df_af.reset_index(drop=True,inplace=True)
df_gt=pd.read_csv(gt,header=None)
df_gt.loc[np.max(df_af, axis=1)>float(cutoff),].to_csv("{}_genotype_matrix.csv".format(cutoff),index=False,header=False)
