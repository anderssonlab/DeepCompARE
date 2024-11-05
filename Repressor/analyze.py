import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu,pearsonr
from loguru import logger


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_two_columns, calc_or, plot_or, plot_or_jitter
from utils import get_track_num



#-------------------Functions-------------------

def read_file(file_suffix):
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv')
    df=df[df['chip_evidence']==True].reset_index(drop=True)
    cols_remove=[col for col in df.columns if col.startswith('pred_orig')]
    df.drop(cols_remove, axis=1, inplace=True)
    df["motif_length"]=df["end"]-df["start"]
    df["max_af"] = df["af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["percentage_variants"] = (df["af"].str.count(":") + 1)/df["motif_length"]
    df["have_common_variant"]= (df["max_af"]>0.001)
    df["241way_max"] = df["241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["447way_max"] = df["447way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["447way_min"] = df["447way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    return df

def add_tf_codependency(df):
    tfs_codependent=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent_merged.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant_merged.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df


def preprocess_group_by_tf(file_suffix):
    df=read_file(file_suffix)
    df=add_tf_codependency(df)
    df=df.groupby("protein").agg({"percentage_variants":"mean", 
                                "have_common_variant":"mean",
                                "241way_max":"mean",
                                "447way_max":"mean",
                                "241way_min":"mean",
                                "447way_min":"mean",
                                }).reset_index()
    return df





#-------------------effect size v.s. max allele frequency -------------------

for file in ["promoters_hepg2", "promoters_k562", "enhancers_hepg2", "enhancers_k562"]:
        df=read_file(file)
        for track_num in get_track_num(file,classification=True):
        # df=df[df["ism_track{track_num}"]<0].reset_index(drop=True)
                df_plot=bin_two_columns(df,
                        f"ism_track{track_num}",
                        [-np.inf,-0.5,-0.2,-0.1, 0, 0.1, 0.2, 0.5, np.inf],
                        "max_af",
                        {"0 - 0.001": "rare", "0.001 - 0.01": "low", "0.01 - 1":"common"},
                        [0,0.001,0.01,1],
                        "variant_type")
                df_plot=df_plot[df_plot.sum(axis=1)>10]
                df_plot=calc_or(df_plot,"Predicted_effect_size","variant_type",out_group="low")
                plot_or(df_plot, 'Predicted_effect_size', 'odds_ratio', "variant_type",
                        f"Odds ratio ({file}, track {track_num})",
                        {'rare': "#1f77b4", 'common': '#ff7f0e'},
                        f"or_max_af_{file}_track{track_num}.pdf")





for file in ["promoters_hepg2", "promoters_k562", "enhancers_hepg2", "enhancers_k562"]:
        df=read_file(file)
        for track_num in get_track_num(file,classification=True):
                df_plot=bin_two_columns(df,
                        f"ism_track{track_num}",
                        [-np.inf, -0.5,-0.2, -0.1, 0, 0.1, 0.2, 0.5, np.inf],
                        "percentage_variants",
                        {"0 - 0.1": "constrained", "0.1 - 0.5": "medium", "0.5 - 2.5":"tolerant"},
                        [0,0.1,0.5,2.5],
                        "TFBS_type")
                df_plot=df_plot[df_plot.sum(axis=1)>10]
                df_plot=calc_or(df_plot,"Predicted_effect_size","TFBS_type",out_group="low")
                plot_or(df_plot, 'Predicted_effect_size', 'odds_ratio', "TFBS_type",
                        f"Odds ratio ({file}, track {track_num})",
                        {'constrained': "#1f77b4", 'tolerant': '#ff7f0e'},
                        f"Plots/or_percentage_variants_{file}_track{track_num}.pdf")




for file in ["promoters_hepg2", "promoters_k562", "enhancers_hepg2", "enhancers_k562"]:
        df=read_file(file)
        #df=df[df["241way_max"]>0].copy()
        for track_num in get_track_num(file,classification=True):
                df_plot=bin_two_columns(df,
                        f"ism_track{track_num}",
                        [-np.inf, -0.5,-0.2, -0.1, 0, 0.1, 0.2, 0.5, np.inf],
                        "241way_max",
                        {"0 - 1": "no_evidence", "1 - 3": "conserved", "3 - inf":"very_conserved"},
                        [0,1,3,np.inf],
                        "TFBS_type")
                df_plot=df_plot[df_plot.sum(axis=1)>10]
                df_plot=calc_or(df_plot,"Predicted_effect_size","TFBS_type",out_group="low")
                plot_or(df_plot, 'Predicted_effect_size', 'odds_ratio', "TFBS_type",
                        f"Odds ratio ({file}, track {track_num})",
                        {'no_evidence': "#1f77b4", 'very_conserved': '#ff7f0e'},
                        f"Plots/or_241way_max_{file}_track{track_num}.pdf")






for file in ["promoters_hepg2", "promoters_k562", "enhancers_hepg2", "enhancers_k562"]:
        df=read_file(file)
        for track_num in get_track_num(file,classification=True):
                df_plot=bin_two_columns(df,
                        f"ism_track{track_num}",
                        [-np.inf, -0.5,-0.2, -0.1, 0, 0.1, 0.2, 0.5, np.inf],
                        "241way_min",
                        {"-1 - 0": "no_evidence", "-3 - -1": "accelerated", "-inf - -3":"very_accelerated"},
                        [-np.inf,-3,-1,0],
                        "TFBS_type")
                df_plot=df_plot[df_plot.sum(axis=1)>10]
                df_plot=calc_or(df_plot,"Predicted_effect_size","TFBS_type",out_group="low")
                plot_or(df_plot, 'Predicted_effect_size', 'odds_ratio', "TFBS_type",
                        f"Odds ratio ({file}, track {track_num})",
                        {'no_evidence': "#1f77b4", 'very_accelerated': '#ff7f0e'},
                        f"Plots/or_241way_min_{file}_track{track_num}.pdf")







