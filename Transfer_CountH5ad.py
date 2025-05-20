import scanpy as sc
import sys

input_file = sys.argv[1]
prefix_input_file = input_file.split(".")[0]

adata = sc.read_h5ad(input_file)
del adata.obs
del adata.var
adata.write("{}_transfer.h5ad".format(prefix_input_file))
