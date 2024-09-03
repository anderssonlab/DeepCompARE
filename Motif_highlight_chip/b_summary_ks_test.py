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
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import match_by_decile

track_info={
    "cage":[0,1],
    "dhs":[2,3],
    "starr":[4,5],
    "sure":[6,7]
            }


#--------------------------------------------------------
# Functions
#--------------------------------------------------------


def ks_with_sign(x,y):
    d_stat,p_val=ks_2samp(x,y)
    return d_stat*np.sign(np.median(x)-np.median(y)), p_val



if __name__ == "__main__":
    
    parser=argparse.ArgumentParser()
    parser.add_argument("--track",type=str,help="Track")
    args=parser.parse_args()
    
    # get empty lists
    tf_list=[]
    occupancy_list=[]
    motif_length_list=[]
    motif_mode_list=[]
    
    sample_size_original_list=[]
    sample_size_balanced_list=[]

    ism_motif_d_stat_list=[]
    ism_motif_p_val_list=[]
    
    mean_ism_motif_list=[]
    median_ism_motif_list=[]
    mean_ism_motif_true_list=[]
    median_ism_motif_true_list=[]

    motif_score_d_stat_list=[]
    motif_score_p_val_list=[]

    # read in data
    df_enhancer_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_enhancers_hepg2_track{track_info[args.track][0]}.csv")
    df_promoters_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_promoters_hepg2_track{track_info[args.track][0]}.csv")
    df_enhancer_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_enhancers_k562_track{track_info[args.track][1]}.csv")
    df_promoters_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_promoters_k562_track{track_info[args.track][1]}.csv")

    # concatinate the 4 dataframes by row, scatter plot the feat_imp_orig vs score, and color by file name (promoters_k562, enhancers_k562, promoters_hepg2)
    df=pd.concat([df_promoters_k562,df_enhancer_k562,df_promoters_hepg2,df_enhancer_hepg2])
    tfs=sorted(df[df.chip_evidence==True].protein.unique())
    logger.info(f"Number of TFs with chip evidence: {len(tfs)}")

    # process each TF
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
        ism_motif_d_stat, ism_motif_p_val=ks_with_sign(df_sub[df_sub.chip_evidence==True].ism_motif,df_sub[df_sub.chip_evidence==False].ism_motif)
        motif_score_d_stat, motif_score_p_val=ks_2samp(df_sub[df_sub.chip_evidence==True].score,df_sub[df_sub.chip_evidence==False].score)
        
        # add to lists
        tf_list.append(tf)
        occupancy_list.append(occupancy)
        
        sample_size_original_list.append(original_sample_size)
        sample_size_balanced_list.append(df_sub.shape[0])
        
        ism_motif_d_stat_list.append(ism_motif_d_stat)
        ism_motif_p_val_list.append(ism_motif_p_val)
        
        mean_ism_motif_list.append(df_sub['ism_motif'].mean())
        median_ism_motif_list.append(df_sub['ism_motif'].median())
        mean_ism_motif_true_list.append(df_sub[df_sub.chip_evidence==True].ism_motif.mean())
        median_ism_motif_true_list.append(df_sub[df_sub.chip_evidence==True].ism_motif.median())
        
        
        motif_score_d_stat_list.append(motif_score_d_stat)
        motif_score_p_val_list.append(motif_score_p_val)

        # violin plot to show distribution difference between chip and non-chip
        sns.violinplot(data=df_sub,x='chip_evidence',y='ism_motif')
        plt.xlabel("Chip evidence")
        plt.ylabel("motif importance")
        plt.title(f"{tf} in {args.track}")
        plt.savefig(os.path.join("Plots_violin_motif_importance",f"{tf}_{args.track}.png"),dpi=300)
        plt.close()
        
    # create report dataframe
    df_ks=pd.DataFrame({"protein":tf_list,
                        "occupancy":occupancy_list,
                        "motif_length":motif_length_list,
                        "motif_mode":motif_mode_list,
                        "ism_motif_d_stat":ism_motif_d_stat_list,
                        "ism_motif_p_val":ism_motif_p_val_list,
                        "mean_ism_motif":mean_ism_motif_list,
                        "median_ism_motif":median_ism_motif_list,
                        "mean_ism_motif_true":mean_ism_motif_true_list,
                        "median_ism_motif_true":median_ism_motif_true_list,
                        "motif_score_d_stat":motif_score_d_stat_list,
                        "motif_score_p_val":motif_score_p_val_list,
                        "sample_size_original":sample_size_original_list,
                        "sample_size_balanced":sample_size_balanced_list})
 
    df_ks.to_csv(f"Pd1_ks_tests/summary_ks_test_{args.track}.csv",index=False)
    
    logger.info(f"Done with {args.track}")
    
    
