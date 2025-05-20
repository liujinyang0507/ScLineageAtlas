import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly as pl
import plotly.express as px
import os
import sys
import pdb

umap_files = sys.argv[1]
clone_files = sys.argv[2]
out_dir = sys.argv[3]

sample = os.path.basename(clone_files).split("_")[0]
print(sample)
df_umap = pd.read_csv(umap_files,header=0,index_col=0)
df_umap["cell type"] = df_umap["cell type"].str.replace(" ","_")
#umap file
#X_umap1    X_umap2         cell_type
#AAACCTGCACACCGAC   2.511441  12.718133        T/NK cells

df_clone = pd.read_csv(clone_files,header=0,index_col=0)
df_clone["clone_id"] = ["clone{}".format(str(i)) for i in df_clone["clone_id"].values]

df = df_umap.merge(df_clone,left_on = df_umap.index,right_on = df_clone.index,how = "inner" )

#df = pd.read_csv("cells_clones_infomation_celltype.csv",index_col=0)
df.loc[:,"celltype_clone"] = df["cell type"].str.cat(df["clone_id"],sep="_")
df.sort_values(by="celltype_clone",axis=0,inplace=True)
print(df.head())

df["celltype_clone"]=df["celltype_clone"].str.replace(" ","_")
df.to_csv(os.path.join(out_dir,"{}_umap_celltype_clone_information.csv".format(sample)),index = False)

0,206,209
colors = [(95,158,160),(0,206,209),(135,206,235),(176,224,230),(255,182,193),(238,130,238),(186,85,211),(147,112,219),(255250,205),(255,160,122),
         (255,192,203),(255,218185),(0,255,255),(250,128,114),(244,164,96),(95,158,160),(48,128,20),(160,32,240),(218,112,214),(240,255,240),(144,238,144),
          (152,251,152),(143,188,143),(50,205,50),(0,255,0),(34,139,34),(0,128,0),(0,100,0),(127,255,0),(124,252,0),(173,255,47),(85,107,47),(245,245,220),(250,250,210)]
cell_types = sorted(df["cell type"].unique())
d_color_map = dict()
for i in range(len(cell_types)):
    df_1 =  df.loc[df["cell type"]== cell_types[i],:]
    clone_ids = sorted(df_1["clone_id"].unique())
    color_rgb = colors[i]
    alpha_values = np.linspace(0.4,1,len(clone_ids)).tolist()
    for j in range(len(clone_ids)):
        df_2 = df_1.loc[df_1["clone_id"]== clone_ids[j]]
        str_celltype_clone = list(set(df_2["celltype_clone"].tolist()))[0]
        color_rgb = list(color_rgb)
        color_rgb.append(alpha_values[j])
        color_rgb_alpha = tuple(color_rgb)
        color_rgb = color_rgb[:3]
        #print({str_celltype_clone:color_rgb_alpha})
        d_color_map[str_celltype_clone] = "rgba{}".format(color_rgb_alpha)
fig = px.scatter(df,x="X_umap1",y="X_umap2",color_discrete_sequence=px.colors.qualitative.Alphabet,
           color = "celltype_clone",
           color_discrete_map = d_color_map,template="simple_white")

pl.offline.plot(fig,filename=os.path.join(out_dir,"{}_umap_clones.html".format(sample)))
plt.axis('off')

df_total=pd.DataFrame(df["cell type"].value_counts())
df_total.columns = ["total"]
clone_id_all = sorted(df["clone_id"].unique())
for i in sorted(clone_id_all):
    df_each = df.loc[df["clone_id"]==i,:]
    df_each = pd.DataFrame(df_each["cell type"].value_counts())
    df_each.columns = [i]
    df_total = pd.concat([df_total,df_each],axis=1)
df_total.fillna(0,inplace=True)

for i in sorted(clone_id_all):
    df_total[i] =df_total.apply(lambda a:a[i]/a["total"],axis=1)
df_total.drop(columns=["total"],inplace=True)

df_total.to_csv(os.path.join(out_dir,"{}_fraction_of_each_clone_of_each_cell_type.csv".format(sample)))

#colours=["cadetblue","sandybrown","salmon","mediumpurple","moccasin","lightgrey","dodgerblue","cadetblue","pink","violet","mediumorchid","skyblue",
#        "cyan","mediumturquoise","aquamarine","mediumspringgreen","palegreen","greenyellow","lemonChiffon","moccasin","tan","peachpuff","lightsalmon"]
colours=['cadetblue', 'darkturquoise','skyblue','powderblue','lightpink','violet', 'mediumorchid','mediumpurple','lemonchiffon','lightsalmon',
         'pink','peachpuff','cyan','salmon','sandybrown','greenyellow','mediumturquoise','lightgrey','palegreen','dodgerblue','mediumspringgreen','tan','aquamarine','aquamarine']
clones = df_total.columns
plt.figure(figsize = (8, 8),dpi=150)
for i in range(len(clones)):
    if i == 0:
        bottom = ""
        plt.bar(x = df_total.index,height = df_total[clones[i]],label = clones[i],color = colours[i],width = 0.5,edgecolor="white")
    elif i == 1:
        bottom = df_total[clones[0]]
        plt.bar(x = df_total.index,height = df_total[clones[i]],bottom = bottom,label = clones[i],color = colours[i],width = 0.5,edgecolor="white")
    else:
        bottom += df_total[clones[i-1]] 
        plt.bar(x = df_total.index,height = df_total[clones[i]],bottom = bottom,label = clones[i],color = colours[i],width = 0.5,edgecolor="white")
plt.xlabel('Cell type',fontsize = 12)   # 横轴标签
plt.ylabel('Fraction',fontsize = 12)  # 纵轴标签
plt.xticks(np.arange(len(cell_types)),cell_types,fontsize = 10,rotation = 45,ha='right')  # 柱状图横轴坐标各类别标签
plt.legend(bbox_to_anchor = (1.03,0.8),title = "clone id",frameon=False)  # 显示两组柱状图的标签
plt.tight_layout()
plt.savefig(os.path.join(out_dir,"{}_fraction_of_each_clone_of_each_cell_type.png".format(sample)))
