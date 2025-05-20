import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
import pandas as pd
import os
from vireoSNP import BinomMixtureVB
import argparse
import seaborn as sns


def FindingOptimalClonesNumber(AD,DP,set_growth_ratio_value,maximum_number_of_clones,min_iter,n_init,out_dir,run_accession):

    AD = mmread(AD).tocsc()
    DP = mmread(DP).tocsc()

    n_clone_list = np.arange(2,int(maximum_number_of_clones) + 2)
    _ELBO_mat_dict =dict()
    _ELBO_mat = []
    df = pd.DataFrame(columns=["n_clone","ELBO_MAX","growth_ratio"])
    df.to_csv(os.path.join(out_dir,
                           "{}_the_evidence_lower_bound_ELBO_values_and_the_number_of_clone_information.csv".format(run_accession)),index = False)
            
    for k in n_clone_list:
        print("set n_clone:{},analysis satrt".format(k))
        _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=k)
        _model.fit(AD, DP, min_iter=min_iter, n_init=n_init)

        _ELBO_mat_dict[k] = max(_model.ELBO_inits)
        _ELBO_mat.append(_model.ELBO_inits)

        if k == 2:
            df = pd.DataFrame({"n_clone":k,"ELBO_MAX":_ELBO_mat_dict[k],"growth_ratio":0},index=[k-1])
            df.to_csv(os.path.join(out_dir,"{}_the_evidence_lower_bound_ELBO_values_and_the_number_of_clone_information.csv".format(run_accession))
                      ,index = False,mode="a",header=False)
            continue
        else:
            growth_ratio = (_ELBO_mat_dict[k] - _ELBO_mat_dict[k-1]) / abs(_ELBO_mat_dict[k-1])
            df = pd.DataFrame({"n_clone":k,"ELBO_MAX":_ELBO_mat_dict[k],"growth_ratio":growth_ratio},index=[k-1])
            df.to_csv(os.path.join(out_dir,"{}_the_evidence_lower_bound_ELBO_values_and_the_number_of_clone_information.csv".format(run_accession)),
                                   index = False,mode="a",header=False)
            if growth_ratio < set_growth_ratio_value:
                n_donor = k - 1
                return n_donor,_ELBO_mat,n_clone_list
            if k == int(int(maximum_number_of_clones) + 1):
                n_donor = k
                return n_donor,_ELBO_mat,n_clone_list

########################################################################
# MAIN
########################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Choose the number of donors based on the evidence lower bound(ELBO) values.'
                                                 'One empirical suggestion is to choose the n_clones when ELBO stops '
                                                 'increasing dramatically')
    parser.add_argument('-ad', '--ad_matrix',
                        help='AD matrix output by MQuad.')
    parser.add_argument('-dp', '--dp_matrix',
                        help='DP matrix output by MQuad.')
    parser.add_argument('-gr', '--growth_ratio_cutoff', type=float,
                        help='The cutoff value of ELBO growth rate',default=0.04)
    parser.add_argument('-mc', '--maximum_number_of_clones', default=20,
                        help='The maximum number of clones')
    parser.add_argument('-o', '--out_dir',
                        help='The result output directory')
    parser.add_argument('-mi', '--min_iter',default=30,type=int,
                        help='The value of min iter')
    parser.add_argument('-ni', '--n_init',default=300,type=int,
                        help='The number of init')
    parser.add_argument('-r', '--run_accession',
                        help='Run Accession')
    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    print("ad_matrix_dir:{}".format(args.ad_matrix))
    print("dp_matrix_dir:{}".format(args.dp_matrix))
    print("The cutoff value of ELBO growth rate:{}".format(args.growth_ratio_cutoff))
    print("The maximum number of clones:{}".format(args.maximum_number_of_clones))
    print("min_iter:{}".format(args.min_iter))
    print("n_init:{}".format(args.n_init))
    print("out_dir:{}".format(args.out_dir))
    print("Run Accession:{}".format(args.run_accession))

    n_donor,_ELBO_mat,n_clone_list = FindingOptimalClonesNumber(args.ad_matrix, args.dp_matrix, args.growth_ratio_cutoff,
                                                                args.maximum_number_of_clones,args.min_iter,args.n_init,
                                                                args.out_dir,args.run_accession)

    pd.DataFrame({"optimal_number_of_clones":n_donor},index=[0]).to_csv(os.path.join(
        args.out_dir,"{}_the_optimal_number_of_clones.csv".format(args.run_accession)),index = False)
    print("The optimal number of clones:{}".format(n_donor))

    df = pd.DataFrame(_ELBO_mat,index=[i for i in range(len(_ELBO_mat))])
    df = df.T
    df.to_csv(os.path.join(args.out_dir, "{}_ELBO.csv".format(args.run_accession)), index=False)
    plt.figure(figsize=(8, 8), dpi=150)
    plt.plot(np.arange(1, n_donor + 1), np.max(df.values, axis=0))
    part = plt.violinplot(df)
    plt.xticks(np.arange(1, n_donor + 1), np.arange(2, n_donor + 2))
    for pc in part["bodies"]:
        print(pc)
        pc.set_facecolor("cadetBlue")
        pc.set_alpha(0.8)
        pc.set_linestyle("--")
    plt.ylabel("ELBO")
    plt.xlabel("n_clones")
    plt.savefig(os.path.join(args.out_dir, "{}_finding_the_optimal_number_of_clones_using_the_ELBO_curve.png".format(args.run_accession)))
    plt.show()

