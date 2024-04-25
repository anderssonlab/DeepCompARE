#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import bin_and_label

#------------------------------------
# Analysis1: plot distribution of feat_imp_orig, split by low/high score
#------------------------------------
df_enhancer_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_no_thresh_enhancers_k562.csv")
df_enhancer_k562["file_name"]="enhancers_k562"

df_enhancer_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_no_thresh_enhancers_hepg2.csv")
df_enhancer_hepg2["file_name"]="enhancers_hepg2"

df_promoters_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_no_thresh_promoters_k562.csv")
df_promoters_k562["file_name"]="promoters_k562"

df_promoters_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_no_thresh_promoters_hepg2.csv")
df_promoters_hepg2["file_name"]="promoters_hepg2"


# concatinate the 3 dataframes by row, scatter plot the feat_imp_orig vs score, and color by file name (promoters_k562, enhancers_k562, promoters_hepg2)
df=pd.concat([df_promoters_k562,df_enhancer_k562,df_promoters_hepg2,df_enhancer_hepg2])
corr, pval = pearsonr(df.feat_imp_orig, df.score)


df=bin_and_label(df, "score", [100,500,1000])

plt.figure(figsize=(12, 6))
sns.violinplot(data=df, x="file_name", y="feat_imp_orig", hue="Bin", inner='quart',split=True)
plt.savefig("Plots/feat_imp_distribution_stratified_by_score.pdf")
plt.xlabel("Jaspar motif score")
plt.ylabel("DeepCompare-assigned motif importance")
plt.close()




#------------------------------------
# Analysis2: calculate correlation between feat_imp_org and "ism_motif"
#------------------------------------
pearsonr(df.feat_imp_orig, df.ism_motif) # 0.61, 0
