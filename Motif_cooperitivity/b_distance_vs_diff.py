""" Understand how TF cooperitivity die down as distance increases"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from loguru import logger

def change_track_name(df):
    track_dict={"track0":"cage","track1":"cage",
                "track2":"dhs","track3":"dhs",
                "track4":"starr","track5":"starr",
                "track6":"sure","track7":"sure"}
    for col in df.columns:
        col_suffix=col.split("_")[-1]
        if col_suffix in track_dict:
            df=df.rename(columns={col:col.replace(col_suffix,track_dict[col_suffix])})
    return df



def rearrange(df):
    df_cage=df[["protein1","protein2","distance","diff_cage"]].copy()
    df_dhs=df[["protein1","protein2","distance","diff_dhs"]].copy()
    df_starr=df[["protein1","protein2","distance","diff_starr"]].copy()
    df_sure=df[["protein1","protein2","distance","diff_sure"]].copy()
    df_cage["data"]="cage"
    df_dhs["data"]="dhs"
    df_starr["data"]="starr"
    df_sure["data"]="sure"
    df_cage=df_cage.rename(columns={"diff_cage":"diff"})
    df_dhs=df_dhs.rename(columns={"diff_dhs":"diff"})
    df_starr=df_starr.rename(columns={"diff_starr":"diff"})
    df_sure=df_sure.rename(columns={"diff_sure":"diff"})
    df=pd.concat([df_cage,df_dhs,df_starr,df_sure])
    return df


def read_file(file_name):
    df=pd.read_csv(file_name)
    if "null" in file_name:
        df=df.drop(columns=["start_rel1","end_rel1"])
        # rename columns
        df=df.rename(columns={"start_rel_null": "start_rel1",
                              "end_rel_null": "end_rel1"})
        for col in df.columns:
            if "null_" in col:
                df=df.rename(columns={col:col.replace("null_","")})
    df["distance"]=np.abs(df["start_rel2"]-df["start_rel1"])
    if "hepg2" in file_name:
        track_list=[0,2,4,6]
    elif "k562" in file_name:
        track_list=[1,3,5,7]
    else:
        raise ValueError("Invalid dataset")
    cols_retain=[f"ism1_track{track_num}" for track_num in track_list]+ \
        [f"ism2_track{track_num}" for track_num in track_list]+ \
        [f"ism_both_track{track_num}" for track_num in track_list]+ \
        ["distance"]
    df=df[["protein1","protein2"]+cols_retain]
    for track_num in track_list:
        df[f"ism2_wo_1_track{track_num}"]=df[f"ism_both_track{track_num}"]-df[f"ism1_track{track_num}"]
        # TODO: is np.abs() needed?
        df[f"diff_track{track_num}"]=df[f"ism2_track{track_num}"]-df[f"ism2_wo_1_track{track_num}"]
    df=change_track_name(df)
    cols_retain=[col for col in df.columns if col.startswith("diff")]
    df=df[["protein1","protein2","distance"]+cols_retain]
    # df=rearrange(df)
    df["avg_diff"]=df[cols_retain].mean(axis=1)
    return df


df_null_enhancers_hepg2=read_file("df_mutate_pair_null_enhancers_hepg2.csv")
df_null_enhancers_k562=read_file("df_mutate_pair_null_enhancers_k562.csv")
df_null_promoters_hepg2=read_file("df_mutate_pair_null_promoters_hepg2.csv")
df_null_promoters_k562=read_file("df_mutate_pair_null_promoters_k562.csv")

df_enhancers_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_hepg2.csv")
df_enhancers_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_k562.csv")
df_promoters_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_hepg2.csv")
df_promoters_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_k562.csv")


df_null=pd.concat([df_null_enhancers_hepg2,df_null_enhancers_k562,df_null_promoters_hepg2,df_null_promoters_k562])
df_null=df_null.groupby(["distance","data"]).agg({"diff":"mean"}).reset_index()

df_alt=pd.concat([df_enhancers_hepg2,df_enhancers_k562,df_promoters_hepg2,df_promoters_k562])
df_alt=df_alt.groupby(["distance","data"]).agg({"diff":"mean"}).reset_index()


sns.lineplot(x="distance", y="diff", data=df_alt, hue="data")
sns.lineplot(x="distance", y="diff", data=df_null, hue="data", alpha=0.3,legend=False)
plt.title("Effect of partner decreases as distance increases")
plt.savefig(f"Plots/distance_vs_diff.pdf")
plt.close()







#---------------------------------------------------
# do co-dependent pairs require closer distance?
#---------------------------------------------------


def ks_test(df):
    df_codependent=df[df["codependent"]]
    df_redundant=df[df["redundant"]]
    d,p=ks_2samp(df_codependent["distance"],df_redundant["distance"])
    sign=df_redundant["distance"].median()-df_codependent["distance"].median()
    d=d*sign
    return d,p,df_redundant["distance"].median(),df_codependent["distance"].median()



df_enhancers_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_hepg2.csv")
df_enhancers_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_k562.csv")
df_promoters_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_hepg2.csv")
df_promoters_k562=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_k562.csv")

df=pd.concat([df_enhancers_hepg2,df_enhancers_k562,df_promoters_hepg2,df_promoters_k562])


# lists to save info
tf_list=[]
dstat_list=[]
pval_list=[]
dataset_list=[]
distance_median_redundant_list=[]
distance_median_codependent_list=[]

dataset_dict={"enhancers_hepg2":df_enhancers_hepg2,
              "enhancers_k562":df_enhancers_k562,
              "promoters_hepg2":df_promoters_hepg2,
              "promoters_k562":df_promoters_k562
              }


for dataset, df in dataset_dict.items():
    df["codependent"]=(df["avg_diff"]>0.01)
    df["redundant"]=(df["avg_diff"]<-0.01)
    for protein in df["protein2"].unique():
        df_sub=df[df["protein2"]==protein].reset_index(drop=True)
        if df_sub.shape[0]<10:
            continue
        if df_sub["codependent"].nunique()==1:
            continue
        # KS test for distance
        dstat,pval,distance_median_redundant,distance_median_codependent=ks_test(df_sub)
        dstat_list.append(dstat)
        pval_list.append(pval)
        distance_median_redundant_list.append(distance_median_redundant)
        distance_median_codependent_list.append(distance_median_codependent)
        print(distance_median_redundant,distance_median_codependent)
        tf_list.append(protein)
        # record the name of the variable
        dataset_list.append(dataset)


df_res=pd.DataFrame({"tf":tf_list,"dstat":dstat_list,"pval":pval_list,
                     "distance_median_redundant":distance_median_redundant_list,
                     "distance_median_codependent":distance_median_codependent_list,
                     "dataset":dataset_list})

df_res.to_csv("df_distance_ks.csv",index=False)
df_res=pd.read_csv("df_distance_ks.csv")
df_res[(df_res["dstat"]>0) & (df_res["pval"]<0.05)].shape[0]





df_res_codependent=df_res[['tf','distance_median_codependent', 'dataset']].copy()
df_res_redundant=df_res[['tf','distance_median_redundant', 'dataset']].copy()
df_res_codependent.rename(columns={"distance_median_codependent":"distance_median"},inplace=True)
df_res_redundant.rename(columns={"distance_median_redundant":"distance_median"},inplace=True)
df_res_codependent["cooperativity"]="codependent"
df_res_redundant["cooperativity"]="redundant"
df_res=pd.concat([df_res_codependent,df_res_redundant])

sns.violinplot(x="dataset",y="distance_median",hue="cooperativity",data=df_res,split=True,inner="quart")
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.28)
plt.title(f"Cooperitivity vs distance")
plt.savefig(f'Plots/distance_median_distribution.pdf')
plt.close()






# nohup python3 distance_vs_diff.py > distance_vs_diff.log &


