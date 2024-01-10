import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy.stats import ks_2samp
from loguru import logger

df=pd.read_csv("motif_df.csv")
tfs=df[df.chip_evidence==True].protein.unique()
logger.info(f"Number of TFs with chip evidence: {len(tfs)}")

#--------------------------------------------------------
# Analysis 2: for each TF, calculate KS D statistics 
# of max_gradxinp distribution grouped by chip_evidence
#--------------------------------------------------------
feat_imp_d_stat_list=[]
feat_imp_p_val_list=[]
feat_imp_true_median_list=[]
feat_imp_false_median_list=[]

motif_score_d_stat_list=[]
motif_score_p_val_list=[]
motif_score_true_median_list=[]
motif_score_false_median_list=[]

i=0
for tf in tfs:
    if i%10==0:
        logger.info(f"{i} tfs done")
        i+=1
    df_sub=df[df.protein==tf].copy().reset_index(drop=True)
    if len(df_sub)==0:
        continue
    feat_imp_d_stat, feat_imp_p_val=ks_2samp(df_sub[df_sub.chip_evidence==True].max_gradxinp,df_sub[df_sub.chip_evidence==False].max_gradxinp)
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
                    "feat_imp_d_stat":feat_imp_d_stat_list,
                    "feat_impp_val":feat_imp_p_val_list,
                    "feat_imp_true_median":feat_imp_true_median_list,
                    "feat_imp_false_median":feat_imp_false_median_list,
                    "motif_score_d_stat":motif_score_d_stat_list,
                    "motif_score_p_val":motif_score_p_val_list,
                    "motif_score_true_median":motif_score_true_median_list,
                    "motif_score_false_median":motif_score_false_median_list
                    })
df_ks.to_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/ks_stat_featimp_and_motif_score.csv",index=False)






#--------------------------------------------------------
# Analysis 1: For each TF, plot the distribution of max_gradxinp, grouped by chip_evidence
#--------------------------------------------------------
# i=0
# for tf in tfs:
#     if i%10==0:
#         logger.info(f"{i} tfs done")
#         i+=1
        
#     df_sub=df[df.protein==tf].copy().reset_index(drop=True)
#     if len(df_sub)==0:
#         continue
#     sns.violinplot(data=df_sub,x='chip_evidence',y='max_gradxinp')
#     plt.xlabel("Chip evidence")
#     plt.ylabel("max_gradxinp")
#     plt.title(f"Jaspar binding sites for {tf}")
#     plt.savefig(os.path.join("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/Plots_violin_max_gradxinp",f"{tf}_violin_max_gradxinp.png"),dpi=300)
#     plt.close()

