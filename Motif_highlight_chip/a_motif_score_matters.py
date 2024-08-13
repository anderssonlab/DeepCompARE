#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, ks_2samp
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label


track_info={"cage":[0,1],"dhs":[2,3],"starr":[4,5],"sure":[6,7]}

def report_kstest(df):
    dstats, pval=ks_2samp(df[df["Bin"]=="100 - 500"].ism_motif, df[df["Bin"]=="500 - 900"].ism_motif)
    logger.info(f"KS test between motif importance (100-500) and motif importance (500-900): {dstats}, {pval}")
    dstats, pval=ks_2samp(df[df["Bin"]=="500 - 900"].ism_motif, df[df["Bin"]=="900 - 1000"].ism_motif)
    logger.info(f"KS test between motif importance (500-900) and motif importance (900-1000): {dstats}, {pval}")
    
    chip_percentage_low=df[df["Bin"]=="100 - 500"].chip_evidence.mean() 
    chip_percentage_medium=df[df["Bin"]=="500 - 900"].chip_evidence.mean() 
    chip_percentage_high=df[df["Bin"]=="900 - 1000"].chip_evidence.mean()
    logger.info(f"Percentage of motifs with ChIP evidence: {chip_percentage_low}, {chip_percentage_medium}, {chip_percentage_high}")
    return


for track in track_info.keys():
    logger.info(f"track: {track}")
    df_enhancer_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_enhancers_hepg2_track{track_info[track][0]}.csv")
    df_enhancer_hepg2["file_name"]="enhancers_hepg2"
    df_promoters_hepg2=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_promoters_hepg2_track{track_info[track][0]}.csv")
    df_promoters_hepg2["file_name"]="promoters_hepg2"

    df_enhancer_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_enhancers_k562_track{track_info[track][1]}.csv")
    df_enhancer_k562["file_name"]="enhancers_k562"
    df_promoters_k562=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_no_thresh_promoters_k562_track{track_info[track][1]}.csv")
    df_promoters_k562["file_name"]="promoters_k562"


    # concatinate the 4 dataframes by row, scatter plot the feat_imp_orig vs score, and color by file name (promoters_k562, enhancers_k562, promoters_hepg2)
    df=pd.concat([df_promoters_k562,df_enhancer_k562,df_promoters_hepg2,df_enhancer_hepg2])
    df["re"]=df["file_name"].apply(lambda x: x.split("_")[0])
    df=bin_and_label(df, "score", [100,500,900,1000])
    corr, pval = pearsonr(df.ism_motif.abs(), df.score) 
    logger.info(f"Correlation between motif importance and Jaspar score: {corr}, p-value: {pval}")
    corr=pearsonr(df.feat_imp_orig, df.ism_motif) 
    logger.info(f"Correlation between gradxinp and ism: {corr}")
    
    report_kstest(df)
    # # subsample df to get equal number of bins
    # df_joint=df.groupby("Bin").apply(lambda x: x.sample(n=min(len(x), 10000), random_state=1))
    # # plot cumulative distribution of motif importance, hue by bin
    # plt.figure(figsize=(6,4))
    # sns.kdeplot(data=df_joint, x="ism_motif", hue="Bin", common_norm=False, cumulative=True)
    # plt.xlabel("Motif importance (DeepCompare)")
    # # change legend title to "Jaspar motif score"
    # plt.savefig(f"Plots/distribution_ism_motif_vs_motif_score_{track}.pdf")
    # plt.close()

    df_promoters=df[df["re"]=="promoters"].reset_index(drop=True)
    logger.info(f"Promoters: {len(df_promoters)}")
    report_kstest(df_promoters)
    # df_promoters=df_promoters.groupby("Bin").apply(lambda x: x.sample(n=min(len(x), 10000), random_state=1))
    df_enhancers=df[df["re"]=="enhancers"].reset_index(drop=True)
    # df_enhancers=df_enhancers.groupby("Bin").apply(lambda x: x.sample(n=min(len(x), 10000), random_state=1))
    logger.info(f"Enhancers: {len(df_enhancers)}")
    report_kstest(df_enhancers)
    # plt.figure(figsize=(6,4))
    # sns.kdeplot(data=df_enhancers, x="ism_motif", hue="Bin", common_norm=False, cumulative=True, linestyle=":")
    # sns.kdeplot(data=df_promoters, x="ism_motif", hue="Bin", common_norm=False, cumulative=True)
    # plt.xlabel("Motif importance (DeepCompare)")
    # plt.title(f"{track}")
    # plt.savefig(f"Plots/distribution_ism_motif_vs_motif_score_{track}_promoters_vs_enhancers.pdf")
    # plt.close()
    logger.info(f"track: {track} done")
    
# python3 a_motif_score_matters.py > a_motif_score_matters.log 2>&1 & 