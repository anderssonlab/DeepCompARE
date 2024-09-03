import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,spearmanr,ks_2samp

import sys

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_file_remove_confusion,detect_and_remove_confusing_pairs



def preprocess_file(file_path,track_nums,remove_confusion=True):
    df=read_file_remove_confusion(file_path,0.1,track_nums=track_nums)
    if remove_confusion:
        df=detect_and_remove_confusing_pairs(df)
        df["codependency"]=(df["codependency_count"]>0).astype(int)
        df.drop(columns=['redundancy_count', 'codependency_count'],inplace=True)
    # rename track0/1 to cage, track2/3 to dhs, track4/5 to starr, track6/7 to sure
    df.rename(columns={"ism2_track0":"cage",
                       "ism2_track1":"cage",
                       "ism2_track2":"dhs",
                       "ism2_track3":"dhs",
                       "ism2_track4":"starr",
                       "ism2_track5":"starr",
                       "ism2_track6":"sure",
                       "ism2_track7":"sure"},inplace=True)
    cols_to_remove=[col for col in df.columns if col.startswith("ci") or col.startswith("cooperativity") or col.startswith("ism")]
    df.drop(columns=cols_to_remove,inplace=True)
    # remove duplicated rows
    df.drop_duplicates(subset=['region_idx','chromosome2','start2','end2','protein2','score2','strand2','start_rel2','end_rel2'],inplace=True)
    return df


df_coop=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_cooperativity_ratio_lenient.csv")

df_promoter_hepg2=preprocess_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_hepg2.csv",track_nums=[0,2,4,6])
df_promoter_k562=preprocess_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_k562.csv",track_nums=[1,3,5,7])
df_enhancer_hepg2=preprocess_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_hepg2.csv",track_nums=[0,2,4,6])
df_enhancer_k562=preprocess_file("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_k562.csv",track_nums=[1,3,5,7])
df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562])
df.fillna(0,inplace=True)


#-------------------------------------------------------------------
# Analysis 1: Correlation between TF effect on different tracks and cooperativity ratio
#-------------------------------------------------------------------

df_grouped=df.groupby(["protein2"]).agg({"cage":"mean","dhs":"mean","starr":"mean","sure":"mean"}).reset_index()
df_grouped.fillna(0,inplace=True)
df_coop=pd.merge(df_grouped,df_coop[["protein2","cooperativity_ratio"]],left_on="protein2",right_on="protein2",how="inner")



# scatter plot
for track in ["cage","dhs","starr","sure"]:
    # small dots
    sns.regplot(x=f"{track}",y="cooperativity_ratio",data=df_coop,scatter_kws={"s":1})
    # annotate with protein name
    for i in range(df_coop.shape[0]):
        plt.text(df_coop[track][i],df_coop["cooperativity_ratio"][i],df_coop["protein2"][i],fontsize=6)
    # annotate with pearson r
    r,p=pearsonr(df_coop[track],df_coop["cooperativity_ratio"])
    plt.text(0.7,0.9,f"pearson r={r:.2f}\np={p:.2f}",transform=plt.gca().transAxes)
    plt.title(f"TF effect on {track} vs cooperativity ratio")
    plt.savefig(f"cr_vs_{track}_lenient.pdf")
    plt.close()






#---------------------------------------------
# Analysis 2: Compare TF effect between redundant and codependent state
#---------------------------------------------


df["avg_effect"]=df[["cage","dhs","starr","sure"]].mean(axis=1)




tf_list=[]
for tf in df["protein2"].unique():
    df_sub=df[df["protein2"]==tf].reset_index(drop=True)
    df_redundant=df_sub[df_sub["codependency"]==0].reset_index(drop=True)
    df_codependent=df_sub[df_sub["codependency"]==1].reset_index(drop=True)
    if df_redundant.shape[0]==0 or df_codependent.shape[0]==0:
        continue
    for track in ["cage","dhs","starr","sure","avg_effect"]:
        stat,p=ks_2samp(df_redundant[track],df_codependent[track])
        # reverse sign if the mean of codependent is smaller
        if df_codependent[track].median()<df_redundant[track].median():
            stat=-stat
        tf_list.append([tf,track,stat,p])

df_result=pd.DataFrame(tf_list,columns=["protein2","track","stat","p"])
# pivot stat
df_stat=df_result.pivot(index="protein2",columns="track",values="stat")
df_pval=df_result.pivot(index="protein2",columns="track",values="p")
df_mask=df_pval>0.05
# apply mask to df_stat
df_stat=df_stat.where(~df_mask.values)
# remove rows with all nan
df_stat.dropna(how="all",inplace=True)
# replace nan with 0
df_stat.fillna(0,inplace=True)

(df_stat["avg_effect"]>0).sum()
(df_stat["cage"]>0).sum()
(df_stat["dhs"]>0).sum()
(df_stat["starr"]>0).sum()
(df_stat["sure"]>0).sum()


df_median=df.groupby(["protein2","codependency"]).agg({"avg_effect":"median"}).reset_index()
df_median=df_median.pivot(index="protein2",columns="codependency",values="avg_effect")
df_median.columns=["redundant_median","codependent_median"]

# merge df_stat with df_median
df_stat=pd.merge(df_stat,df_median,left_on="protein2",right_on="protein2",how="inner")
df_stat["significant"]=(df_stat["avg_effect"]!=0)

# merge with cooperativity ratio
df_stat=pd.merge(df_stat,df_coop[["protein2","cooperativity_ratio"]],left_on="protein2",right_on="protein2",how="inner")


# any rows with all negative values?
df_stat[(df_stat[["cage","dhs","starr","sure"]]<0).all(axis=1)]
df_stat[(df_stat[["cage","dhs","starr","sure"]]>=0).all(axis=1)]


# shape by significance, hue by cooperativity ratio
sns.scatterplot(x="redundant_median",y="codependent_median",data=df_stat,hue="cooperativity_ratio",style="significant")
plt.title("tf effect comparison")
# add diagonal line
min_val=min(df_stat["redundant_median"].min(),df_stat["codependent_median"].min())
max_val=max(df_stat["redundant_median"].max(),df_stat["codependent_median"].max())
plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
plt.savefig("redundant_median_vs_codependent_median_lenient.png")
plt.close()








