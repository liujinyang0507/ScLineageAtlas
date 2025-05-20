import re
import pandas as pd
import sys

input_file = sys.argv[1]
out_file = sys.argv[2]

df=pd.read_csv(input_file,header=0)
def replace_str(str):
    res=re.search(r"(.*?)([AGCT]{8,}(-1)?)",str)
    child_str1 = res.group(1)
    child_str2 = res.group(2)
    if child_str1:
        new_str = re.sub("-", "_", child_str1) + child_str2
        return new_str
    else:
        new_str = child_str2
        return new_str
df["new_index"] = df.iloc[:,0].apply(replace_str)
df.drop(df.columns[0], axis=1,inplace=True)
df_ = df.loc[:,["new_index","clone_id"]]
df_.columns=["barcodes","clone_id"]
df_.to_csv(out_file,index=False)