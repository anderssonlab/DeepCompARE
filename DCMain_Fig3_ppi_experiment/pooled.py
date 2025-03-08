import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer
from tf_cooperativity import assign_cooperativity


ci_suffix="pe"
redundancy_threshold=0.3 # 0.46
codependent_threshold=0.7 # 0.83
mode="linearity_index"

df_coop=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/Temp/tf_pair_cooperativity_index_k562_{ci_suffix}.csv")
if mode=="cooperativity_index":
    df_coop=df_coop[df_coop["c_sum"]>1].reset_index(drop=True)
if mode=="linearity_index":
    df_coop=assign_cooperativity(df_coop,redundancy_threshold,codependent_threshold)
df_coop=df_coop[df_coop["protein2"].isin(["BACH1","MAFG","IKZF1","RREB1","RFX5"])].reset_index(drop=True)



# pool
df_bait_pooled=pd.DataFrame()
for bait in ["BACH1","MAFG","IKZF1","RREB1"]:
    df_bait=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
    df_bait["bait"]=bait
    df_bait_pooled=pd.concat([df_bait_pooled,df_bait],axis=0)


# rename df_bait_pooled
df_bait_pooled.rename(columns={"gene":"protein1","bait":"protein2"},inplace=True)
# merge
df_bait_pooled=pd.merge(df_bait_pooled,df_coop,how="inner",on=["protein1","protein2"])


df_failed_pooled=pd.DataFrame()
for bait in ["BACH1","MAFG","IKZF1","RREB1","RFX5"]:
    df_failed=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.FilteredProteins.txt", sep='\t')
    df_failed["reference"]=df_failed["reference"].apply(lambda x: x.split("|")[-1])
    df_failed["gene"]=df_failed["reference"].apply(lambda x: x.split("_")[0])
    df_failed["species"]=df_failed["reference"].apply(lambda x: x.split("_")[-1])
    # remove nonhuman
    df_failed=df_failed[df_failed["species"]=="HUMAN"].reset_index(drop=True)
    df_failed["bait"]=bait
    df_failed_pooled=pd.concat([df_failed_pooled,df_failed],axis=0)


df_failed_pooled.rename(columns={"gene":"protein1","bait":"protein2"},inplace=True)
df_failed_pooled=pd.merge(df_failed_pooled,df_coop,how="inner",on=["protein1","protein2"])




# plot distribution of cooperativity index
sns.kdeplot(df_coop[mode],color="blue",cut=0,label="All")
sns.kdeplot(df_bait_pooled[mode],color="red",cut=0,label="TFs reported")
sns.kdeplot(df_failed_pooled[mode],color="green",cut=0,label="TFs filtered")
plt.xlabel(mode)
plt.ylabel("Density")
plt.title(f"Distribution of {mode} ({ci_suffix})")
plt.legend()
plt.savefig(f"pooled_{mode}_distribution_{ci_suffix}.pdf")
plt.close()

