# coding: utf-8

from scipy.io import mmread
from scipy.io import mmwrite
from scipy import sparse
import pandas as pd
import argparse
import pdb

def split_AD_DP_files_to_each_celltype_by_barcode(AD_matrix,DP_matrix,dataset_barcodes,cell_types,out_dir):
    AD = mmread(AD_matrix).tocsc()
    DP = mmread(DP_matrix).tocsc()

    df_AD = pd.DataFrame(AD.todense())
    df_DP = pd.DataFrame(DP.todense())


    df_dataset_barcodes = pd.read_csv(dataset_barcodes,header=None)
    df_cell_types = pd.read_csv(cell_types,header=None)
    cell_types_list = df_cell_types.iloc[:,0].tolist()
    print(cell_types_list)
    for cell_type in cell_types_list:
        #pdb.set_trace()
        cell_type_index = df_dataset_barcodes.loc[df_dataset_barcodes.iloc[:,1] == cell_type,:].index.tolist()

        df_AD_ = df_AD.iloc[:,cell_type_index]
        df_AD_.reset_index(drop = True,inplace = True)
        #AD_sparse = sparse.coo_matrix((df_AD.values))
        mmwrite(out_dir + '/{}_passed_ad.mtx'.format(cell_type), sparse.csr_matrix(df_AD_))

        df_DP_ = df_DP.iloc[:,cell_type_index]
        df_DP_.reset_index(drop = True,inplace = True)
        #DP_sparse = sparse.coo_matrix((df_DP.values))
        mmwrite(out_dir + '/{}_passed_dp.mtx'.format(cell_type), sparse.csr_matrix(df_DP_))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='split AD and DP files to each celltype by barcodes')
    parser.add_argument('-ad', '--AD',
                        help='A sparse matrix of AD in a mtx format')
    parser.add_argument('-dp', '--DP',
                        help='A sparse matrix of DP in a mtx format')
    parser.add_argument('-bd', '--dataset_barcodes',
                        help='A list containing barcodes and cell type')
    parser.add_argument('-ct', '--cell_type_list',
                        help='A list 0f cell type')
    parser.add_argument('-o', '--out_dir',
                        help='out dir')
    args = parser.parse_args()

    split_AD_DP_files_to_each_celltype_by_barcode(args.AD,args.DP,args.dataset_barcodes,args.cell_type_list,args.out_dir)






