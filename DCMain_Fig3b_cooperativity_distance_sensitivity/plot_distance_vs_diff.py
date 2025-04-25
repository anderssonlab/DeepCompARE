""" Understand how TF cooperitivity die down as distance increases"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import change_track_name


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


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
    cols_retain=[f"isa1_track{track_num}" for track_num in track_list]+ \
        [f"isa2_track{track_num}" for track_num in track_list]+ \
        [f"isa_both_track{track_num}" for track_num in track_list]+ \
        ["distance"]
    df=df[["protein1","protein2"]+cols_retain]
    for track_num in track_list:
        df[f"isa2_wo_1_track{track_num}"]=df[f"isa_both_track{track_num}"]-df[f"isa1_track{track_num}"]
        df[f"diff_track{track_num}"]=np.abs(df[f"isa2_track{track_num}"]-df[f"isa2_wo_1_track{track_num}"])
    df=change_track_name(df)
    cols_retain=[col for col in df.columns if col.startswith("diff")]
    df=df[["protein1","protein2","distance"]+cols_retain]
    df=rearrange(df)
    return df






#option=""
option="_cnn6"

df_null_enhancers_hepg2=read_file(f"Pd1_mutate_pair_null/df_mutate_pair_null{option}_enhancers_hepg2.csv")
df_null_enhancers_k562=read_file(f"Pd1_mutate_pair_null/df_mutate_pair_null{option}_enhancers_k562.csv")
df_null_promoters_hepg2=read_file(f"Pd1_mutate_pair_null/df_mutate_pair_null{option}_promoters_hepg2.csv")
df_null_promoters_k562=read_file(f"Pd1_mutate_pair_null/df_mutate_pair_null{option}_promoters_k562.csv")

df_enhancers_hepg2=read_file(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs{option}_enhancers_hepg2.csv")
df_enhancers_k562=read_file(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs{option}_enhancers_k562.csv")
df_promoters_hepg2=read_file(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs{option}_promoters_hepg2.csv")
df_promoters_k562=read_file(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs{option}_promoters_k562.csv")


df_null=pd.concat([df_null_enhancers_hepg2,df_null_enhancers_k562,df_null_promoters_hepg2,df_null_promoters_k562])
df_null=df_null.groupby(["distance","data"]).agg({"diff":"mean"}).reset_index()

df_alt=pd.concat([df_enhancers_hepg2,df_enhancers_k562,df_promoters_hepg2,df_promoters_k562])
df_alt=df_alt.groupby(["distance","data"]).agg({"diff":"mean"}).reset_index()

plt.figure(figsize=(2.3,2))
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)

sns.lineplot(x="distance", y="diff", data=df_alt, hue="data", linewidth=0.5)
sns.lineplot(x="distance", y="diff", data=df_null, hue="data", legend=False, linewidth=0.5, linestyle="--")

plt.xlabel("Distance", fontsize=7)
plt.ylabel("Interaction", fontsize=7)

# Set log scale and specify ticks
plt.xscale("log")
xticks = [10, 20, 50, 100, 200, 300, 600]
plt.xticks(xticks, [str(x) for x in xticks], fontsize=5)

y_ticks = [0, 0.05, 0.1, 0.15, 0.2]
plt.yticks(y_ticks, fontsize=5)

plt.legend(title="Modality", fontsize=5, title_fontsize=5)
plt.tight_layout()
plt.savefig(f"Plots/distance_vs_diff{option}.pdf")
plt.close()




# nohup python3 plot_distance_vs_diff.py > plot_distance_vs_diff.log &


