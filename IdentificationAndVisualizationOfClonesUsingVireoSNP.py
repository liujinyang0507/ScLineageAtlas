import vireoSNP
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
import matplotlib as mpl
from vireoSNP import BinomMixtureVB
import seaborn as sns

#from mquad.mquad_utils import plot_confusionMatrix, confusionMatrix
from vireoSNP.plot.base_plot import heat_matrix
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
import os
import argparse

print(vireoSNP.__version__)
np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})

def IdentificationAndVisualizationOfClones(AD,DP,variants,barcode,n_donor,min_iter,n_init,out_dir,run_accession):
    mquad_AD = mmread(AD).tocsc()
    mquad_DP = mmread(DP).tocsc()
    variants = pd.read_table(variants, header=None)
    cells_num = mquad_AD.shape[1]
    with open(barcode, 'r') as f:
        sample_id = f.read().splitlines()

    #fit on mquad output
    np.random.seed(42)
    _model = BinomMixtureVB(n_var=len(mquad_AD.getnnz(axis=1)), n_cell=len(mquad_AD.getnnz(axis=0)), n_donor=n_donor)
    _model.fit(mquad_AD, mquad_DP, min_iter=min_iter, n_init=n_init)

    mquad_modelCA = _model

    ## how many cells are assignable?
    a = np.sum(np.max(mquad_modelCA.ID_prob, axis=1) > 0.8)/cells_num
    print("Percentage of assignable cells: ", a)

    df_assignment_prob = pd.DataFrame(mquad_modelCA.ID_prob,columns=["clone{}".format(i) for i in range(n_donor)],
                                      index=sample_id)
    df_assignment_prob.sort_values(by=df_assignment_prob.columns.tolist(), inplace=True, ascending=False)
    df_assignment_prob.reset_index(drop=True, inplace=True)
    df_assignment_prob.to_csv(os.path.join(out_dir,"{}_cell_assignment_probability_to_each_donor.csv".format(run_accession)))

    #assignment prob heatmap
    plt.figure(figsize=(8, 8), dpi=150)
    im = heat_matrix(mquad_modelCA.ID_prob, cmap="Oranges", alpha=0.8,
                     display_value=False, row_sort=True, interpolation='none')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title("Assignment probability")
    plt.xlabel("Clone")
    plt.ylabel("Cell counts")
    plt.xticks(range(mquad_modelCA.n_donor))
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,"{}_cell_assignment_probability_to_each_donor.png".format(run_accession)))
    plt.savefig(os.path.join(out_dir,"{}_cell_assignment_probability_to_each_donor.pdf".format(run_accession)))

    #Mean allelic ratio heatmap
    top = cm.get_cmap('Blues', 200)
    newcolors = np.vstack((top(np.linspace(0, 0.7, 10)),
                           top(np.linspace(0.7, 1, 90))))
    newcmp = ListedColormap(newcolors, name='segBlues')
    AF_SNPs = mquad_modelCA.beta_mu
    # rearrange clones to match fig 3b
    # AF_SNPs_sorted = np.array([[i[1], i[2], i[0]]for i in AF_SNPs])
    df_mean_allelic_ratio = pd.DataFrame(AF_SNPs,columns=[i for i in range(n_donor)],index = variants.iloc[:,0].tolist())
    df_mean_allelic_ratio.to_csv(os.path.join(out_dir,"{}_mean_allelic_ratio.csv".format(run_accession)))
    plt.figure(figsize=(8, 8), dpi=150)
    sns.heatmap(df_mean_allelic_ratio, cmap=sns.light_palette("cadetblue", as_cmap=True))
    plt.title("Mean allelic ratio")
    plt.xlabel("Clone")
    plt.ylabel("Clonal mtSNV variants")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir,"{}_mean_allelic_ratio.png".format(run_accession)))
    plt.savefig(os.path.join(out_dir,"{}_mean_allelic_ratio.pdf".format(run_accession)))

    #plot the AF heatmap to visualize their differences between clones
    AF_df = pd.DataFrame(mquad_AD/mquad_DP, index = variants.iloc[:,0].tolist(),columns=sample_id).fillna(0)
    AF_df.to_csv(os.path.join(out_dir,"{}_allele_frequency.csv".format(run_accession)))

    ## perform hclust on rows instead of using snsclustermap cuz its slow af

    from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
    linked = linkage(AF_df, 'single')
    plt.subplot(1,1,1)
    dendrogram(linked,
                orientation='top',
                distance_sort='descending',
                labels=AF_df.index,
                show_leaf_counts=True)
    plt.xticks(rotation=45, fontsize=4)
    plt.savefig(os.path.join(out_dir,"Allele_frequency_cluster.svg"),format='svg')
    plt.savefig(os.path.join(out_dir,"Allele_frequency_cluster.pdf"),format='pdf')

    cmap_c = plt.get_cmap("tab20c")
    cmap_b = plt.get_cmap("tab20b")
    new_cmap = ListedColormap(cmap_c.colors+cmap_b.colors)

    ## AF matrix
    top = cm.get_cmap('Greens', 200)

    newcolors = np.vstack((top(np.linspace(0, 0.7, 10)),
                           top(np.linspace(0.7, 1, 90))))
    newcmp = ListedColormap(newcolors, name='segGreens')

    plt.figure(figsize=(8, 8), dpi=150)
    ax = plt.subplot(111)
    clone_id = np.argmax(mquad_modelCA.ID_prob, axis=1)

    clones_df = pd.DataFrame(data={'sample_id': sample_id, 'clone_id': clone_id})
    clones_df.to_csv(os.path.join(out_dir, "{}_cells_clones_infomation.csv".format(run_accession)),index=False)

    row_idx = leaves_list(linked)
    res = ax.imshow(AF_df.iloc[row_idx, np.argsort(clone_id)], cmap=newcmp, aspect='auto', interpolation='none')
    ax.axes.xaxis.set_visible(False)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    plt.yticks(range(len(AF_df)), variants.iloc[row_idx][0], fontsize=7)
    divider = make_axes_locatable(ax)
    left_ax = divider.append_axes("left", size=0.5, pad=0)
    top_ax = divider.append_axes("top", size=0.13, pad=0.02, sharex=ax)

    #ax_t = fig.add_axes([0.026,0.07,0.07,0.55], frame_on=False)
    with plt.rc_context({'lines.linewidth': 0.5}):
        d = dendrogram(linked, ax=left_ax, orientation='left', no_labels=True, color_threshold=0, link_color_func=lambda x: 'k')
    left_ax.invert_yaxis()
    left_ax.axis('off')

    #top_ax = fig.add_axes([0.095, 0.65, 0.38, 0.03], frame_on=False)
    #dont change label num
    sorted_df = AF_df.iloc[row_idx, np.argsort(clone_id)].T.reset_index(drop=True).T

    i = 0
    dict = {}
    LEFT=0
    label_num=np.bincount(clone_id)
    for num in label_num:
        print("i, num, left = ", i, num, LEFT)
        dict[str(i)] = sorted_df.loc[:,LEFT:LEFT+num].mean(axis=1)
        print(dict[str(i)].shape)
        LEFT += num
        i += 1

    grouped = pd.DataFrame.from_dict(dict)
    genotype = grouped.apply(lambda x: pd.cut(x, bins=[0, 0.01, 0.9, 1, 5], labels=['0', '1', '2', '3'], include_lowest=True))
    genotype.to_csv(os.path.join(out_dir,'{}_genotype_matrix.csv'.format(run_accession)), sep= ' ', header=False, index=False)

    LEFT = 0
    iicolor = 0
    c_names = [i for i in range(n_donor)]
    label_num=np.bincount(clone_id)
    #clone_pal = ['royalblue', 'darkgreen', 'firebrick', 'purple', 'yellow']
    clone_pal = new_cmap.colors[:n_donor]
    for num in label_num:
            top_ax.barh(0,num+1,left=LEFT,color=clone_pal[iicolor])
            top_ax.text(x=LEFT + num/2, y=0.8, s=c_names[iicolor], va='center', ha='center')
            top_ax.set_xlim(0,cells_num)
            top_ax.axis('off')
            LEFT += num
            iicolor += 1

    bottom_ax = divider.append_axes("bottom", size=0.1, pad = 0.3)
    plt.colorbar(res, cax = bottom_ax, orientation="horizontal", shrink=0.5)
    plt.savefig(os.path.join(out_dir,"{}_allele_frequency_heatmap.png".format(run_accession)))
    plt.savefig(os.path.join(out_dir,"{}_allele_frequency_heatmap.pdf".format(run_accession)))

# MAIN
########################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identification and visualization of clones using vireoSNP\
                                                 by manually setting the number of clones')
    parser.add_argument('-ad', '--ad_matrix',
                        help='AD matrix output by MQuad.')
    parser.add_argument('-dp', '--dp_matrix',
                        help='DP matrix output by MQuad.')
    parser.add_argument('-var', '--passed_varients',
                        help='passed varients list output by MQuad.')
    parser.add_argument('-b', '--barcode',
                        help='barcode files')
    parser.add_argument('-n', '--donor_counts',type=int,
                        help='Manually setting the number of clones')
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
    print("passed_varients_dir:{}".format(args.passed_varients))
    print("barcode_dir:{}".format(args.barcode))
    print("The number of clones was set:{}".format(args.donor_counts))
    print("min_iter:{}".format(args.min_iter))
    print("n_init:{}".format(args.n_init))
    print("out_dir:{}".format(args.out_dir))
    print("Run Accession:{}".format(args.run_accession))

    IdentificationAndVisualizationOfClones(args.ad_matrix, args.dp_matrix, args.passed_varients,args.barcode,
                                           args.donor_counts,args.min_iter,args.n_init,args.out_dir,args.run_accession)





