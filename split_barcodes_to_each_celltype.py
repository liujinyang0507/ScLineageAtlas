import pandas as pd
import os
import sys
import re
import argparse
import pdb

# Match the corresponding samples based on the SRR number
def split_obs(obs_file,sample_ino_file,patient,out_file):
    df_obs = pd.read_csv(obs_file,header=0,index_col=0)
    df_barcode = pd.read_csv(sample_ino_file,header=0)

    # try:
    sample = df_sample_ino.loc[df_sample_ino["run_accession"]==patient,"sample_title"].tolist()[0]
    barcodes_list = list()
    for i in df_obs.index.tolist():
        #print(i)
        barcode = re.search("([AGCT]{12,}-*\d*)", i).group(1)
        text = re.compile(r".*[0-9]$")
        if text.match(barcode):
            barcodes_list.append(barcode)
        else:
            barcodes_list.append("{}-1".format(barcode))

    df_obs["barcodes"] = barcodes_list

    colnames_obs = df_obs.columns.tolist()
    for col in colnames_obs:
        res = re.search("(sample.*)", col,re.I)
        if res:
            df_obs = df_obs.loc[:, [res.group(1), "cell type", "barcodes"]]
            df_obs_1 = df_obs[df_obs[res.group(1)] == sample]
            break
        else:
            res = re.search("(patient.*)", col, re.I)
            if res:
                df_obs = df_obs.loc[:, [res.group(1), "cell type", "barcodes"]]
                df_obs_1 = df_obs[df_obs[res.group(1)] == sample]
                break

    df_obs_1.loc[:,"cell type"] =  df_obs_1["cell type"].str.replace(" ","_")
    df_obs_1.to_csv(out_file,header=None)
    # except Exception as e:
    #     print("Error,because of {}".format(e))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split cell for each sample')
    parser.add_argument('-b', '--obs',
                        help='A file containing barcodes,paper patient ID and cell type')
    parser.add_argument('-e', '--sample_ino_file',
                        help='A file containing paper patient ID and SRR ID,coluns=["run_accession","sample_title"]')
    parser.add_argument('-s', '--sample_id',
                        help='sample SRR id')
    parser.add_argument('-o', '--out_file',
                        help='out file')
    args = parser.parse_args()
    print("obs_file:{}".format(args.obs))
    print("sample_ino_file:{}".format(args.sample_ino_file))
    print("sample_id:{}".format(args.sample_id))
    print("out_file:{}".format(args.out_file))
    split_obs(args.obs, args.sample_ino_file, args.sample_id, args.out_file)


import os
os.chdir("")
df_obs = pd.read_csv(obs_file,header=0,index_col=0)
df_barcode = pd.read_csv(sample_ino_file,header=0)

# try:
sample = df_sample_ino.loc[df_sample_ino["run_accession"]==patient,"sample_title"].tolist()[0]
barcodes_list = list()
for i in df_obs.index.tolist():
    #print(i)
    barcode = re.search("([AGCT]{12,}-*\d*)", i).group(1)
    text = re.compile(r".*[0-9]$")
    if text.match(barcode):
        barcodes_list.append(barcode)
    else:
        barcodes_list.append("{}-1".format(barcode))

df_obs["barcodes"] = barcodes_list

colnames_obs = df_obs.columns.tolist()
for col in colnames_obs:
    res = re.search("(sample.*)", col,re.I)
    if res:
        df_obs = df_obs.loc[:, [res.group(1), "cell type", "barcodes"]]
        df_obs_1 = df_obs[df_obs[res.group(1)] == sample]
        break
    else:
        res = re.search("(patient.*)", col, re.I)
        if res:
            df_obs = df_obs.loc[:, [res.group(1), "cell type", "barcodes"]]
            df_obs_1 = df_obs[df_obs[res.group(1)] == sample]
            break

df_obs_1.loc[:,"cell type"] =  df_obs_1["cell type"].str.replace(" ","_")
df_obs_1.to_csv(out_file,header=None)
# except Exception as e:
#     print("Error,because of {}".format(e))





