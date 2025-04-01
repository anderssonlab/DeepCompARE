import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

from scipy.stats import pearsonr
from scipy.stats import fisher_exact



re="distal"


def group_mediator_interactors_per_region(df_med,mediator,df_tfbs):
    if mediator is not "all":
        df_med=df_med[df_med["bait"]==mediator].reset_index(drop=True)
    # is protein in df_med.gene
    df=df_tfbs.copy()
    df[mediator]=df["protein"].apply(lambda x: x in df_med["gene"].values.tolist())
    # group by region, and count the number of mediator, get the mean gc content
    df=df.groupby("region").agg({mediator:"sum"}).reset_index()
    df[mediator] = pd.cut(df[mediator], bins=[-1, 0, 1, 2, 3, 100], labels=["0", "1", "2", "3", "3+"])
    return df



def calc_percentage_above_thresh_signal(df,mediator,col):
    # get 75 percentile
    thresh_signal=df["log1p_"+col].quantile(0.75)
    df["above_thresh_signal"]=df["log1p_"+col]>thresh_signal
    # fisher exact of above median signal within each mediator-interactor group
    contingency_table=pd.crosstab(df[mediator],df["above_thresh_signal"])
    # add column "percentage>median"
    contingency_table["percentage>median"]=contingency_table[True]/(contingency_table[True]+contingency_table[False])
    contingency_table["mediator"]=mediator
    return contingency_table




def get_ctables(df,df_med,col):
    ctables=pd.DataFrame()
    for mediator in df_med["bait"].unique().tolist()+["all"]:
        df_temp=calc_percentage_above_thresh_signal(df,mediator,col)
        if ctables.empty:
            ctables=df_temp
        else:
            ctables=pd.concat([ctables,df_temp],axis=0)
    return ctables.reset_index()




def plot(ctables,col):
    plt.figure(figsize=(4, 3))
    sns.lineplot(x="index",y="percentage>median",data=ctables,marker="o",hue="mediator")
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel("# mediator-interactor TFBSs",fontsize=7)
    plt.ylabel("fraction of above median regions",fontsize=7)
    plt.title(f"{col},{re}",fontsize=7)
    plt.tight_layout()
    plt.legend(fontsize=5)
    plt.savefig(f"{col}_{re}.pdf")
    plt.close()




# read histone marks
df_histone=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/{re}_k562.csv")
# TODO: choose to -1 or not
df_histone["end"]=df_histone["end"]
# add region in format chr:start-end
df_histone["region"]=df_histone["chrom"].astype(str)+":"+df_histone["start"].astype(str)+"-"+df_histone["end"].astype(str)



# # plot correlation between histone marks
# corr_map=df_histone.iloc[:,3:].corr()
# # rename columns: remove log1p_
# corr_map.columns=[col.replace("log1p_","") for col in corr_map.columns]
# corr_map.index=[col.replace("log1p_","") for col in corr_map.index]
# # plot clustermap
# plt.figure()
# sns.clustermap(corr_map)
# plt.savefig(f"clustermap_histone_{re}.pdf")
# plt.close()





# read tfbs
df_tfbs=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs_{re}_k562.csv")
# subset for proper background TFs
proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
df_tfbs=df_tfbs[df_tfbs['protein'].isin(proteins_background)].reset_index(drop=True)


# read mediator interactors
df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)

# output the "gene" in text format
import numpy as np
pd.Series(np.sort(df_med["gene"].unique())).to_csv("mediator_interactors.txt",index=False,header=False)


# count number of mediator-interactor TFBSs per region
df_counts=pd.DataFrame()
for mediator in df_med["bait"].unique().tolist()+["all"]:
    df_temp = group_mediator_interactors_per_region(df_med, mediator, df_tfbs)
    # if df_counts is empty, set df_counts to df_res
    if df_counts.empty:
        df_counts=df_temp
    else:
        df_counts=pd.merge(df_temp,df_counts,on="region",how="inner")

# merge df_counts with df_histone
df=pd.merge(df_counts,df_histone,on="region",how="inner")






log1p_cols = [col for col in df_histone.columns if "log_signal" in col]
log1p_cols = [col.replace("log1p_", "") for col in log1p_cols]



for col in log1p_cols:
    ctables=get_ctables(df,df_med,col)
    plot(ctables,col)



























#------------------------------------------------
#  avg number of TFBSs per promoter v.s. enhancer
#------------------------------------------------

# df_promoter=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_k562.csv")
# df_enhancer=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_k562.csv")

# df_med=pd.read_csv(f"2024-08-07_MED-TF_interactions.txt",sep="\t")
# df_med=df_med[df_med["significant"]].reset_index(drop=True)
# tfs=df_med["gene"].unique().tolist()

# # select only tfs
# df_promoter=df_promoter[df_promoter["protein"].isin(tfs)].reset_index(drop=True)
# df_enhancer=df_enhancer[df_enhancer["protein"].isin(tfs)].reset_index(drop=True)

# tf_counts_promoter=df_promoter.protein.value_counts()
# tf_counts_enhancer=df_enhancer.protein.value_counts()

# # how many promoter regions in total?
# n_promoter_regions=len(df_promoter["region"].unique())
# n_enhancer_regions=len(df_enhancer["region"].unique())

# # normalize by number of regions
# tf_counts_promoter=tf_counts_promoter/n_promoter_regions
# tf_counts_enhancer=tf_counts_enhancer/n_enhancer_regions

# # merge tf counts by name 
# df=pd.merge(tf_counts_promoter,tf_counts_enhancer,left_index=True,right_index=True,how="inner",suffixes=("_promoter","_enhancer"))

# # rename protein_promoter to promoter
# df.rename(columns={"protein_promoter":"promoter","protein_enhancer":"enhancer"},inplace=True)

# # scatter plot 
# plt.figure()
# plt.scatter(df["promoter"],df["enhancer"])
# # add name for each point
# for i in range(len(df)):
#     plt.text(df["promoter"][i],df["enhancer"][i],df.index[i])

# # add diagonal line
# min_val=min(df["promoter"].min(),df["enhancer"].min())
# max_val=max(df["promoter"].max(),df["enhancer"].max())
# plt.plot([min_val,max_val],[min_val,max_val],color="gray",linestyle="--",linewidth=0.5)
# plt.xlabel("avg # TFBSs per promoter region")
# plt.ylabel("avg # TFBSs per enhancer region")
# plt.title("avg # TFBSs per region")
# plt.savefig("avg_tfbs_per_region.png")
# plt.close()














# nohup python3 tfbs.py > tfbs.out &