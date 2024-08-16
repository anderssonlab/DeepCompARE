import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_promoter_vs_enhancer(df,cell_type,ism):
    df_sub=df[df["dataset"].str.contains(cell_type)].reset_index(drop=True)
    df_sub=df_sub[['protein',f'dstat_ism_{ism}', 'dataset']]
    df_sub=df_sub.pivot(index='protein', columns='dataset', values=f'dstat_ism_{ism}').reset_index()
    df_sub["opposite_sign"]=df_sub[f"promoters_{cell_type}"]*df_sub[f"enhancers_{cell_type}"]<0
    df_sub=df_sub[df_sub["opposite_sign"]]
    df_promoter=df_sub[df_sub[f"promoters_{cell_type}"]-df_sub[f"enhancers_{cell_type}"]>0.2]
    df_enhancer=df_sub[df_sub[f"enhancers_{cell_type}"]-df_sub[f"promoters_{cell_type}"]>0.2]
    tfs_promoter=df_promoter["protein"].tolist()
    tfs_enhancer=df_enhancer["protein"].tolist()
    return tfs_promoter,tfs_enhancer


def dict2df(this_dict):
    unique_tfs = sorted(list(set(sum(this_dict.values(), []))))
    df = pd.DataFrame(index=unique_tfs, columns=this_dict.keys())
    for col in this_dict:
        df[col] = df.index.isin(this_dict[col]).astype(int)
    return df


def promoter_enhancer_contribution_to_tf_specificity(cell_type):
    df=pd.read_csv("tf_individual_effect_by_file.csv")
    dict_promoter_highlight =  {"cage": [], "dhs": [], "starr": [], "sure": []}
    dict_enhancer_highlight = {"cage": [], "dhs": [], "starr": [], "sure": []}
    for ism in ["cage","dhs","starr","sure"]:
        tfs_promoter,tfs_enhancer=get_promoter_vs_enhancer(df,cell_type,ism)
        dict_promoter_highlight[ism] = tfs_promoter
        dict_enhancer_highlight[ism] = tfs_enhancer
    df_promoter_highlight = dict2df(dict_promoter_highlight)
    df_enhancer_highlight = dict2df(dict_enhancer_highlight)
    df_cell_type_specific=pd.read_csv(f"tfs_{cell_type}_highlight.csv",index_col=0)
    df_cell_type_specific["num_support"]=df_cell_type_specific.sum(axis=1)
    df_cell_type_specific=df_cell_type_specific[df_cell_type_specific["num_support"]>=2]
    set_promoter=set(df_promoter_highlight.index)
    set_enhancer=set(df_enhancer_highlight.index)
    set_both=set(df_cell_type_specific.index)
    print("Promoter",len(set_promoter.intersection(set_both))/len(set_both))
    print("Enhancer",len(set_enhancer.intersection(set_both))/len(set_both))


promoter_enhancer_contribution_to_tf_specificity("hepg2")
promoter_enhancer_contribution_to_tf_specificity("k562")



def remove_both_negatives(df):
    idx_remove=df[(df["enhancers"]<0) & (df["promoters"]<0)].index
    df=df[~df.index.isin(idx_remove)]
    return df


def split_hepg2_k562(df):
    hepg2_cols=[col for col in df.columns if "hepg2" in col]
    df_hepg2=df[hepg2_cols].copy()
    df_hepg2.rename(columns={col:col.replace("_hepg2","") for col in df_hepg2.columns},inplace=True)
    k562_cols=[col for col in df.columns if "k562" in col]
    df_k562=df[k562_cols].copy()
    df_k562.rename(columns={col:col.replace("_k562","") for col in df_k562.columns},inplace=True)
    df_hepg2=remove_both_negatives(df_hepg2)
    df_k562=remove_both_negatives(df_k562)
    return df_hepg2,df_k562



df=pd.read_csv("tf_individual_effect_by_file.csv")
usfs=pd.read_csv("universal_stripe_factors.txt",sep="\t",header=None)
# subset to contain only universal stripe factors
usfs.columns=["protein"]
df=df[df["protein"].isin(usfs["protein"])]

def plot_usfs_enhancer_vs_promoter(df,ism):
    df_sub=df[["protein","dataset",f"dstat_ism_{ism}"]]
    df_sub=df_sub.pivot(index='protein', columns='dataset', values=f'dstat_ism_{ism}')
    df_hepg2,df_k562=split_hepg2_k562(df_sub)
    plt.figure(figsize=(5,5))
    sns.scatterplot(data=df_hepg2,x="enhancers",y="promoters",label="HepG2")
    sns.scatterplot(data=df_k562,x="enhancers",y="promoters",label="K562")
    plt.xlabel("Importance in Enhancer")
    plt.ylabel("Importance in Promoter")
    plt.legend()
    # add abline
    x_min=min(df_hepg2["enhancers"].min(),df_k562["enhancers"].min())
    x_max=max(df_hepg2["enhancers"].max(),df_k562["enhancers"].max())
    y_min=min(df_hepg2["promoters"].min(),df_k562["promoters"].min())
    y_max=max(df_hepg2["promoters"].max(),df_k562["promoters"].max())
    min_val=min(x_min,y_min)
    max_val=max(x_max,y_max)
    plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
    # add horizontal and vertical line
    plt.axhline(0,color="black",linestyle="--")
    plt.axvline(0,color="black",linestyle="--")
    # annotate every dot
    for i, txt in enumerate(df_hepg2.index):
        plt.annotate(txt, (df_hepg2["enhancers"].iloc[i], df_hepg2["promoters"].iloc[i]),fontsize=6)
    for i, txt in enumerate(df_k562.index):
        plt.annotate(txt, (df_k562["enhancers"].iloc[i], df_k562["promoters"].iloc[i]),fontsize=6)
    plt.savefig(f"Plots/usfs_enhancers_vs_promoters_{ism}.pdf",bbox_inches='tight')
    plt.close()



plot_usfs_enhancer_vs_promoter(df,"cage")
plot_usfs_enhancer_vs_promoter(df,"dhs")
plot_usfs_enhancer_vs_promoter(df,"starr")
plot_usfs_enhancer_vs_promoter(df,"sure")











#-----------
# Archived
#-----------
# def plot_dict(this_dict,title,output_file):
#     df=dict2df(this_dict)
#     df = df.reindex(sorted(df.columns), axis=1)
#     plt.figure(figsize=(3, 20))
#     my_plot = sns.heatmap(df, cmap="vlag", linewidths=0.5, annot=False)
#     my_plot .set_xticklabels(my_plot .get_xticklabels(), rotation=90)
#     plt.title(title)
#     plt.xticks(rotation=90)
#     plt.tight_layout()
#     plt.savefig(output_file, dpi=300, bbox_inches='tight')
#     plt.close()


# plot_dict(dict_promoter_highlight,"TFs highlighted by promoter","Plots/tfs_promoters_highlight.pdf")
# plot_dict(dict_enhancer_highlight, "TFs highlighted by enhancer","Plots/tfs_enhancer_highlight.pdf")
# plot_dict(dict_both_highlight, "TFs highlighted in both promoters and enhancers","Plots/tfs_both_re_highlight.pdf")

