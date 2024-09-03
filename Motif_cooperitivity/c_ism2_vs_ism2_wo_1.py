import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


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



def read_file(file_name):
    df=pd.read_csv(file_name)
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
        df[f"diff_track{track_num}"]=np.abs(df[f"ism2_track{track_num}"]-df[f"ism2_wo_1_track{track_num}"])
    df=change_track_name(df)
    return df



df_enhancers_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_hepg2.csv")
df_promoters_hepg2=read_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_hepg2.csv")
df_all=pd.concat([df_enhancers_hepg2,df_promoters_hepg2],axis=0)


tf="HNF4A"
suffix="sure"
df=df_promoters_hepg2
df_sub=df[df["protein2"]==tf].copy()
plt.figure(figsize=(5,5))
# smaller scatter plot
sns.scatterplot(data=df_sub, x=f'ism2_wo_1_{suffix}', y=f'ism2_{suffix}',hue="distance",s=5)
min_val=min(df_sub[f'ism2_wo_1_{suffix}'].min(),df_sub[f'ism2_{suffix}'].min())
max_val=max(df_sub[f'ism2_wo_1_{suffix}'].max(),df_sub[f'ism2_{suffix}'].max())
df_codependent=df_sub[df_sub[f'ism2_wo_1_{suffix}']<df_sub[f'ism2_{suffix}']]
df_redundant=df_sub[df_sub[f'ism2_wo_1_{suffix}']>df_sub[f'ism2_{suffix}']]
# annotate number of codependent pairs
plt.text(min_val+0.1,max_val-0.1,f"# codependent pairs: {df_codependent.shape[0]}",fontsize=8)
plt.text(min_val+0.2,min_val+0.1,f"# redundant pairs: {df_redundant.shape[0]}",fontsize=8)
plt.plot([min_val,max_val],[min_val,max_val],color='black',linestyle='--',linewidth=0.5)
plt.title(f"{tf}, promoters")
plt.legend(loc='lower right')
plt.savefig(f"Plots/{tf}_track_{suffix}_promoters.pdf",dpi=300)
plt.close()





