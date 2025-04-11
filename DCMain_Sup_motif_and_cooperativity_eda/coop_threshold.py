import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

def read_file(file_name):
    df=pd.read_csv(file_name)
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
    return df


df=read_file('/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_hepg2.csv')
# histogram on diff

# rearrange
df_track0=df[['protein1', 'protein2','diff_track0']]
df_track0.columns=['protein1', 'protein2','diff']
df_track0['track']='track0'
df_track2=df[['protein1', 'protein2','diff_track2']]
df_track2.columns=['protein1', 'protein2','diff']
df_track2['track']='track2'
df_track4=df[['protein1', 'protein2','diff_track4']]
df_track4.columns=['protein1', 'protein2','diff']
df_track4['track']='track4'
df_track6=df[['protein1', 'protein2','diff_track6']]
df_track6.columns=['protein1', 'protein2','diff']
df_track6['track']='track6'
df_diff=pd.concat([df_track0,df_track2,df_track4,df_track6])


sns.histplot(data=df_diff, x="diff", hue="track", element="step", common_norm=False,alpha=0.5)
# y log scale
plt.yscale('log')
plt.title("Promoters HepG2")
plt.savefig("diff_promoters_hepg2.png")
plt.close()
