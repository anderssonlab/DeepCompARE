import pandas as pd
import numpy as np
import os
import csv
from loguru import logger
from communities.algorithms import louvain_method


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_annotators import find_binding_evidence



def read_file_remove_confusion(file_name,tpm_thresh, nonlinearity_thresh,track_nums=None):
    """
    Read in file, remove mistaking rows, calculate cooperativity, and return the dataframe with cooperativity information
    """
    # read file
    logger.info(f"Reading {file_name}")
    df=pd.read_csv(file_name)
    logger.info(f"# Total tf pairs: {df.shape[0]}")
    # remove rows if protein1 or protein2 are not expressed
    if "hepg2" in file_name:
        expressed_protein_list=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2_{tpm_thresh}.tsv", sep='\t', header=None).iloc[:,0].tolist()
    elif "k562" in file_name:
        expressed_protein_list=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562_{tpm_thresh}.tsv", sep='\t', header=None).iloc[:,0].tolist()

    df=find_binding_evidence(df,"protein1",expressed_protein_list,"rna_evidence")
    df=df[df["rna_evidence"]].reset_index(drop=True)
    df=find_binding_evidence(df,"protein2",expressed_protein_list,"rna_evidence")
    df=df[df["rna_evidence"]].reset_index(drop=True)
    
    logger.info(f"# Expressed tf pairs: {df.shape[0]}")
    
    tracks_to_remove=[i for i in range(16) if i not in track_nums]
    # remove all columns with track{i} and i is in tracks_to_remove
    col_suffix_to_remove=[f"_track{i}" for i in tracks_to_remove]
    col_names_to_remove=[col for col in df.columns if any([col.endswith(suffix) for suffix in col_suffix_to_remove])]
    df.drop(columns=col_names_to_remove,inplace=True)
    # calculate cooperativity index
    for i in track_nums:
        df[f"isa2_wo_protein1_track{i}"]=df[f'isa_both_track{i}']-df[f'isa1_track{i}']
        # remove potentially mistaking rows and conditional repressors
        df.loc[(df[f"isa1_track{i}"]<0) | (df[f"isa2_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df.loc[(df[f"isa2_wo_protein1_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df[f"c_track{i}"]=df[f"isa2_track{i}"]-df[f"isa2_wo_protein1_track{i}"]
        # assign cooperativity
        df[f"cooperativity_track{i}"]="unknown"
        df.loc[df[f"c_track{i}"]< -nonlinearity_thresh,f"cooperativity_track{i}"]="redundancy"
        df.loc[df[f"c_track{i}"]> nonlinearity_thresh,f"cooperativity_track{i}"]="codependency"
    # drop columns starting with "pred" and isa, and motif_length
    col_names_to_remove=[col for col in df.columns if col.startswith("pred") or col.startswith("motif_length")]
    df.drop(columns=col_names_to_remove,inplace=True)
    return df




def _sum_cooperativity_by_row(df):
    cooperativity_cols=[col for col in df.columns if "cooperativity" in col]
    df["redundancy_count"]=0
    df["codependency_count"]=0
    for col in cooperativity_cols:
        df["redundancy_count"]+= (df[col]=="redundancy")
        df["codependency_count"]+= (df[col]=="codependency")
    return df





def detect_and_remove_confusing_pairs(df):
    df=_sum_cooperativity_by_row(df)
    df=df[(df["redundancy_count"]>0) | (df["codependency_count"]>0)].reset_index(drop=True)
    logger.info(f"# Total tf pairs: {df.shape[0]}")
    row_idx=df[(df["redundancy_count"]>0) & (df["codependency_count"]>0)].index
    logger.info(f"# Confusing tf pairs: {len(row_idx)}")
    df.drop(row_idx,inplace=True)
    df=df.reset_index(drop=True)
    return df




def read_cooperativity(file_name,tpm_thresh, nonlinearity_thresh,track_nums=None):
    df=read_file_remove_confusion(file_name,tpm_thresh,nonlinearity_thresh,track_nums)
    df=detect_and_remove_confusing_pairs(df)
    df.drop(columns=["redundancy_count","codependency_count"],inplace=True)
    df["distance"]=np.abs(df["start1"]-df["start2"])
    # count number of NaN per row
    cooperativity_cols=[col for col in df.columns if "cooperativity" in col]
    df["num_unknown"]=df[cooperativity_cols].apply(lambda row: sum(row=="unknown"),axis=1)
    # which num_nan is NaN
    df["num_valid_profiles"]=len(cooperativity_cols)-df["num_unknown"]
    df.fillna(0,inplace=True)
    # aggregate cooperativity index
    ci_cols=[col for col in df.columns if "c_track" in col]
    df["c"]=df[ci_cols].sum(axis=1)/df["num_valid_profiles"]
    df["codependency"]=(df["c"]>0).astype(int)
    df.drop(columns=ci_cols+cooperativity_cols+["num_unknown","num_valid_profiles"],inplace=True)
    return df



def calculate_tf_pair_cooperativity_index(df):
    """
    df should be the c of various tf pairs
    For each tf pair,
    summarize the c into redundancy and codependency,
    calculate cooperativity_index
    """
    df["distance_iqr"]=df["distance"]
    # count num of rows in each group
    df1=df.groupby(["protein1","protein2","codependency"]).agg({"c":"sum",
                                                               "chromosome1":"count"                                                       
                                                               }).reset_index()
    df1.rename(columns={"chromosome1":"count"},inplace=True)
    df1["tf_pair"]=df1["protein1"]+"_"+df1["protein2"]
    df_pivot=df1.pivot(index="tf_pair",columns="codependency",values=["c","count"]).reset_index()
    df_pivot.columns = ['_'.join(map(str, col)).strip('_') for col in df_pivot.columns]
    df_pivot.columns = [col.replace('_0', '_redundancy').replace('_1', '_codependency') for col in df_pivot.columns]
    df_pivot["protein1"]=df_pivot["tf_pair"].apply(lambda x: x.split("_")[0])
    df_pivot["protein2"]=df_pivot["tf_pair"].apply(lambda x: x.split("_")[1])
    df_pivot.drop("tf_pair",axis=1,inplace=True)
    df_pivot.fillna(0,inplace=True)
    df2=df.groupby(["protein1","protein2"]).agg({"distance":"median",
                                                "distance_iqr": lambda x: np.percentile(x, 75) - np.percentile(x, 25)
                                                }).reset_index()
    # merge df1 and df2
    df_pivot=pd.merge(df_pivot,df2,on=["protein1","protein2"],how="left")
    df_pivot=df_pivot.loc[:,['protein1', 'protein2', 'c_redundancy', 'c_codependency', 'count_redundancy','count_codependency', 'distance','distance_iqr']]
    # calculate sum of cooperativity
    df_pivot["c_sum"]=df_pivot["c_redundancy"].abs()+df_pivot["c_codependency"]
    df_pivot["cooperativity_index"]=df_pivot["c_codependency"].abs()/df_pivot["c_sum"]
    df_pivot["cooperativity_fraction"]=df_pivot["count_codependency"]/(df_pivot["count_redundancy"]+df_pivot["count_codependency"])
    return df_pivot





def calculate_tf_cooperativity_index(df):
    """
    df should be tf pair cooperativity index (pre-filter)
    """
    df_tf=df.groupby("protein2").agg({"c_redundancy":"sum",
                                      "c_codependency":"sum",
                                      "count_redundancy":"sum",
                                      "count_codependency":"sum",
                                      }).reset_index()
    df_tf["c_sum"]=df_tf["c_redundancy"].abs()+df_tf["c_codependency"]
    df_tf["count_sum"]=df_tf["count_redundancy"]+df_tf["count_codependency"]
    df_tf["cooperativity_index"]=df_tf["c_codependency"]/(df_tf["c_redundancy"].abs()+df_tf["c_codependency"])
    df_tf["cooperativity_fraction"]=df_tf["count_codependency"]/df_tf["count_sum"]
    return df_tf



def make_symmetric(df):
    df["pair"]=df.apply(lambda x: "-".join(sorted([x["protein1"],x["protein2"]])),axis=1)
    df_grouped=df.groupby("pair").agg({"cooperativity_index":"mean"}).reset_index()
    df1=df_grouped.copy()
    df2=df_grouped.copy()
    df1[['protein1','protein2']] = df1['pair'].str.split('-',expand=True)
    df2[['protein2','protein1']] = df2['pair'].str.split('-',expand=True)
    df=pd.concat([df1,df2],axis=0)
    df=df.drop_duplicates().reset_index(drop=True)
    df=df.drop(columns=["pair"])
    return df


def get_ppi(df,cooperativity,thresh,file_name):
    """
    df: tf pair cooperativity ratio
    """
    # make df symmetric
    df=make_symmetric(df)
    df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
    protein_names=df.index.tolist()
    if cooperativity=="redundancy":
        mat=(df<thresh).astype(int)
    if cooperativity=="codependency":
        mat=(df>thresh).astype(int)
    community_idx, _ = louvain_method(mat.values)
    res=[]
    for this_set in community_idx:
        this_set=list(this_set)
        res_set=[]
        for i in range(len(this_set)):
            res_set.append(protein_names[this_set[i]])
        if len(res_set)>1:
            res.append(res_set)
    with open(file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for row in res:
            writer.writerow(row)
    return res

