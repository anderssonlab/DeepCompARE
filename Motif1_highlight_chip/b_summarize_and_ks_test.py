"""
This module calculates occupancy and KS statistics for each TF
Using data /isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{}.csv
"""



import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

from scipy.stats import ks_2samp
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from stat_tests import match_by_decile



#--------------------------------------------------------
# Functions
#--------------------------------------------------------

def ks_with_sign(x,y):
    d_stat,p_val=ks_2samp(x,y)
    return d_stat*np.sign(np.median(x)-np.median(y)), p_val



if __name__ == "__main__":
    
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=str,help="Input file prefix")
    args=parser.parse_args()
    
    # get empty lists
    tf_list=[]
    occupancy_list=[]
    motif_length_list=[]
    motif_mode_list=[]
    
    sample_size_original_list=[]
    sample_size_balanced_list=[]

    feat_imp_d_stat_list=[]
    feat_imp_p_val_list=[]
    
    mean_feat_imp_list=[]
    median_feat_imp_list=[]
    mean_feat_imp_true_list=[]
    median_feat_imp_true_list=[]
    
    remove_context_d_stat_list=[]
    remove_context_p_val_list=[]

    motif_score_d_stat_list=[]
    motif_score_p_val_list=[]

    # read in data
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{args.input}.csv")
    tfs=sorted(df[df.chip_evidence==True].protein.unique())
    logger.info(f"Number of TFs with chip evidence: {len(tfs)}")

    # process each TF
    # TODO: add occupancy information for all TFs
    for tf in tfs:
        logger.info(f"Process {tf}")

        df_sub=df[df.protein==tf].copy().reset_index(drop=True)
        if len(df_sub)==0:
            continue
        
        if len(df_sub.chip_evidence.unique())==1:
            continue
        
        occupancy=np.sum(df_sub.chip_evidence)/len(df_sub)
        # adjust score difference
        original_sample_size=df_sub.shape[0]
        df_sub=match_by_decile(df_sub, score_col='score', label_col='chip_evidence')
        
        if df_sub.shape[0]<20:
            continue
        
        # motif length
        df_sub['motif_length']=df_sub['end']-df_sub['start']+1
        motif_length_list.append(df_sub['motif_length'].mean())
        motif_mode_list.append(df_sub['motif_sequence'].mode().values[0])
        
        # ks test on feature importance
        feat_imp_d_stat, feat_imp_p_val=ks_with_sign(df_sub[df_sub.chip_evidence==True].feat_imp_orig,df_sub[df_sub.chip_evidence==False].feat_imp_orig)
        remove_context_d_stat, remove_context_p_val=ks_with_sign(df_sub.feat_imp_orig,df_sub.feat_imp_remove_context)
        motif_score_d_stat, motif_score_p_val=ks_2samp(df_sub[df_sub.chip_evidence==True].score,df_sub[df_sub.chip_evidence==False].score)
        
        # add to lists
        tf_list.append(tf)
        occupancy_list.append(occupancy)
        
        sample_size_original_list.append(original_sample_size)
        sample_size_balanced_list.append(df_sub.shape[0])
        
        feat_imp_d_stat_list.append(feat_imp_d_stat)
        feat_imp_p_val_list.append(feat_imp_p_val)
        
        mean_feat_imp_list.append(df_sub['feat_imp_orig'].mean())
        median_feat_imp_list.append(df_sub['feat_imp_orig'].median())
        mean_feat_imp_true_list.append(df_sub[df_sub.chip_evidence==True].feat_imp_orig.mean())
        median_feat_imp_true_list.append(df_sub[df_sub.chip_evidence==True].feat_imp_orig.median())
        
        remove_context_d_stat_list.append(remove_context_d_stat)
        remove_context_p_val_list.append(remove_context_p_val)
        
        motif_score_d_stat_list.append(motif_score_d_stat)
        motif_score_p_val_list.append(motif_score_p_val)

        # violin plot to show distribution difference between chip and non-chip
        sns.violinplot(data=df_sub,x='chip_evidence',y='feat_imp_orig')
        plt.xlabel("Chip evidence")
        plt.ylabel("motif importance")
        plt.title(f"{tf} in {args.input}")
        plt.savefig(os.path.join("Plots_violin_motif_importance",f"{tf}_{args.input}.png"),dpi=300)
        plt.close()
        
    # create report dataframe
    df_ks=pd.DataFrame({"protein":tf_list,
                        "conditional_occupancy":occupancy_list,
                        "motif_length":motif_length_list,
                        "motif_mode":motif_mode_list,
                        "feat_imp_d_stat":feat_imp_d_stat_list,
                        "feat_imp_p_val":feat_imp_p_val_list,
                        "mean_feat_imp":mean_feat_imp_list,
                        "median_feat_imp":median_feat_imp_list,
                        "mean_feat_imp_true":mean_feat_imp_true_list,
                        "median_feat_imp_true":median_feat_imp_true_list,
                        "remove_context_d_stat":remove_context_d_stat_list,
                        "remove_context_p_val":remove_context_p_val_list,
                        "motif_score_d_stat":motif_score_d_stat_list,
                        "motif_score_p_val":motif_score_p_val_list,
                        "sample_size_original":sample_size_original_list,
                        "sample_size_balanced":sample_size_balanced_list})
 
    df_ks.to_csv(f"summary_and_ks_test_{args.input}.csv",index=False)
    
    logger.info(f"Done with {args.input}")
