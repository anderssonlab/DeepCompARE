#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

from scipy.spatial import cKDTree
from scipy.stats import ks_2samp,pearsonr
from adjustText import adjust_text
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import match_by_score

df=pd.read_csv("motif_df.csv")
tfs=sorted(df[df.chip_evidence==True].protein.unique())
logger.info(f"Number of TFs with chip evidence: {len(tfs)}")


#--------------------------------------------------------
# Analysis 1: 
# for each TF, calculate KS D statistics of max_gradxinp distribution and motif score
# grouped by chip_evidence
#--------------------------------------------------------


corr, p_val = pearsonr(df.max_gradxinp, df.score)   

occupancy_orig_list=[]
occupancy_matched_list=[]

motif_length_list=[]
motif_mode_list=[]

feat_imp_d_stat_list=[]
feat_imp_p_val_list=[]
feat_imp_true_median_list=[]
feat_imp_false_median_list=[]

motif_score_d_stat_list=[]
motif_score_p_val_list=[]
motif_score_true_median_list=[]
motif_score_false_median_list=[]


for tf in tfs:
    logger.info(f"Process {tf}")

    df_sub=df[df.protein==tf].copy().reset_index(drop=True)
    if len(df_sub)==0:
        continue
    
    occupancy_orig_list.append(np.sum(df_sub.chip_evidence)/len(df_sub))
    df_sub=match_by_score(df_sub, score_col='score', label_col='chip_evidence')
    occupancy_matched_list.append(np.sum(df_sub.chip_evidence)/len(df_sub))
    
    # motif length
    df_sub['motif_length']=df_sub['end']-df_sub['start']+1
    motif_length_list.append(df_sub['motif_length'].mean())
    motif_mode_list.append(df_sub['motif_sequence'].mode().values[0])
    
    # plot
    sns.violinplot(data=df_sub,x='chip_evidence',y='max_gradxinp')
    plt.xlabel("Chip evidence")
    plt.ylabel("max_gradxinp")
    plt.title(f"Jaspar binding sites for {tf}")
    plt.savefig(os.path.join("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/Plots_violin_max_gradxinp",f"{tf}_violin_max_gradxinp.png"),dpi=300)
    plt.close()
    
    # compute ks stat
    feat_imp_d_stat, feat_imp_p_val=ks_2samp(df_sub[df_sub.chip_evidence==True].max_gradxinp,df_sub[df_sub.chip_evidence==False].max_gradxinp)
    
    # add to li
    feat_imp_d_stat_list.append(feat_imp_d_stat)
    feat_imp_p_val_list.append(feat_imp_p_val)
    feat_imp_true_median_list.append(df_sub[df_sub.chip_evidence==True].max_gradxinp.median())
    feat_imp_false_median_list.append(df_sub[df_sub.chip_evidence==False].max_gradxinp.median())
    
    motif_score_d_stat, motif_score_p_val=ks_2samp(df_sub[df_sub.chip_evidence==True].score,df_sub[df_sub.chip_evidence==False].score)
    motif_score_d_stat_list.append(motif_score_d_stat)
    motif_score_p_val_list.append(motif_score_p_val)
    motif_score_true_median_list.append(df_sub[df_sub.chip_evidence==True].score.median())
    motif_score_false_median_list.append(df_sub[df_sub.chip_evidence==False].score.median())

df_ks=pd.DataFrame({"protein":tfs,
                    "occupancy_orig":occupancy_orig_list,
                    "occupancy_matched":occupancy_matched_list,
                    "motif_length":motif_length_list,
                    "motif_mode":motif_mode_list,
                    "feat_imp_d_stat":feat_imp_d_stat_list,
                    "feat_impp_val":feat_imp_p_val_list,
                    "feat_imp_true_median":feat_imp_true_median_list,
                    "feat_imp_false_median":feat_imp_false_median_list,
                    "motif_score_d_stat":motif_score_d_stat_list,
                    "motif_score_p_val":motif_score_p_val_list,
                    "motif_score_true_median":motif_score_true_median_list,
                    "motif_score_false_median":motif_score_false_median_list
                    })

# convert all p-values to -log10(p)
df_ks['feat_impp_val']=-np.log10(df_ks['feat_impp_val']+1e-300)
df_ks['motif_score_p_val']=-np.log10(df_ks['motif_score_p_val']+1e-300)

# if true median is less than false median, then the d_stat and -log(p) should be negative
df_ks['feat_imp_d_stat']=df_ks['feat_imp_d_stat']*np.sign(df_ks['feat_imp_true_median']-df_ks['feat_imp_false_median'])
df_ks['feat_impp_val']=df_ks['feat_impp_val']*np.sign(df_ks['feat_imp_true_median']-df_ks['feat_imp_false_median'])
df_ks['motif_score_d_stat']=df_ks['motif_score_d_stat']*np.sign(df_ks['motif_score_true_median']-df_ks['motif_score_false_median'])
df_ks['motif_score_p_val']=df_ks['motif_score_p_val']*np.sign(df_ks['motif_score_true_median']-df_ks['motif_score_false_median'])  

df_ks.to_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/ks_stat_featimp_and_motif_score_v2.csv",index=False)


#--------------------------------------------------------
# Analysis 2:
# Plot occupacy and feat_imp_d_stat
#--------------------------------------------------------

df_ks=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/ks_stat_featimp_and_motif_score_v2.csv")
np.sum(df_ks.motif_score_d_stat>0) # 0, perfect match of scores
np.sum(df_ks.feat_impp_val>-np.log10(0.05)) # 159 TFs with ks test p<0.05
# sort by d_stat
df_ks=df_ks.sort_values(by=['feat_imp_d_stat'],ascending=True).reset_index(drop=True)
corr, p_val = pearsonr(df_ks.feat_imp_d_stat, df_ks.occupancy_orig) # 0.47, 2e-12
corr, p_val = pearsonr(df_ks.feat_imp_d_stat, df_ks.motif_length) # 0.10, 0.14

texts = []
plt.figure(figsize=(10, 8))  
scatter_plot = sns.scatterplot(data=df_ks, x='feat_imp_d_stat', y='occupancy_orig')
for line in range(0, df_ks.shape[0]):
    text=scatter_plot.text(df_ks.feat_imp_d_stat[line], df_ks.occupancy_orig[line], 
                            df_ks.protein[line], horizontalalignment='left', 
                            size='small', color='black')
    texts.append(text)
adjust_text(texts)
plt.title("KS d statistics of feat_imp and occupancy")
plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/Plots_ks/after_score_matching_ks_dstat_featimp_and_occupancy.pdf",dpi=300)
plt.close()



texts = []
plt.figure(figsize=(10, 8))  
scatter_plot = sns.scatterplot(data=df_ks, x='feat_imp_d_stat', y='motif_length')
for line in range(0, df_ks.shape[0]):
    text=scatter_plot.text(df_ks.feat_imp_d_stat[line], df_ks.motif_length[line], 
                            df_ks.protein[line], horizontalalignment='left', 
                            size='small', color='black')
    texts.append(text)
adjust_text(texts)
plt.title("KS d statistics of feat_imp and motif_length")
plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/Plots_ks/after_score_matching_ks_dstat_featimp_and_motif_length.pdf",dpi=300)
plt.close()





## Todo adjust for GC content.
## Work in TF families,
## Case study: why do TFs with same motif bind differently?
