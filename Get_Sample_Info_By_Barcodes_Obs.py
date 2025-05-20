import pandas as pd
import os
import re
from collections import Counter
import sys

sample = sys.argv[1]
obs = sys.argv[2]
barcode = sys.argv[3]
out_file = sys.argv[4]

df_obs=pd.read_csv(obs,header=0,index_col=0)
df_barcode=pd.read_csv(barcode,header=None)

barcodes_raw=df_barcode.iloc[:,0].tolist()

barcodes_list = list()
for i in df_obs.index.tolist():
    barcode = re.search("([AGCT]{12,}-*\d*)", i).group(1)
    text = re.compile(r".*[0-9]$")
    if text.match(barcode):
        barcodes_list.append(barcode)
    else:
        barcodes_list.append("{}-1".format(barcode))
df_obs["barcodes"] = barcodes_list

df_obs__filter_sample = df_obs.loc[df_obs["barcodes"].isin(barcodes_raw),"barcodes"]

new_df = pd.DataFrame(df_obs__filter_sample)

new_df.loc[:,"barcodes_clean"]=new_df.index.tolist()

patients = list()
for i in range(new_df.shape[0]):
    patient = new_df.iloc[i,1].replace("_{}".format(new_df.iloc[i,0]),"")
    patients.append(patient)

# Count the occurrences of each element in the list
counter = Counter(patients)
# Get the element with the most occurrences
patient=counter.most_common(1)[0][0]
if len(set(patients))==1:
    patient2=""
else:
    patient2=counter.most_common(2)[1][0]    
df_sample_patient = pd.DataFrame({"run_accession":sample,"sample_title":patient,"sample_title_2":patient2},index=[0])
print({sample:counter.most_common(len(set(patients)))})
df_sample_patient.to_csv(out_file,header=None,index=False)


