#!/usr/bin/python3
import torch
import numpy as np
from prediction import compute_predictions
from region_ops import resize_df
from loguru import logger
from utils import find_available_gpu


def predict_vcf(df_vcf,seq_extractor,device=torch.device("cuda:"+find_available_gpu())):
    # change default device to predetermined device
    if isinstance(device,str):
        device=torch.device(f"cuda:{device}")
    df_copy=df_vcf.copy()
    df_copy=resize_df(df_copy,600)
    df_copy["seq_ref"]=df_copy.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]-1), axis=1)
    
    # get the center base, validate the center base is same as REF
    df_copy["seq_ref_center"]=df_copy["seq_ref"].apply(lambda x: x[300])
    assert np.all(df_copy["REF"]==df_copy["seq_ref_center"])
    
    # change center base to ALT
    df_copy["prefix"]=df_copy["seq_ref"].apply(lambda x: x[0:300])
    df_copy["suffix"]=df_copy["seq_ref"].apply(lambda x: x[301:600])
    df_copy["seq_alt"]=df_copy["prefix"]+df_copy["ALT"]+df_copy["suffix"]   
    
    # validate seq_alt and seq_ref differ only at center base
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[300])!=df_copy["seq_ref"].apply(lambda x: x[300]))
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[299])==df_copy["seq_ref"].apply(lambda x: x[299]))
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[301])==df_copy["seq_ref"].apply(lambda x: x[301]))
    
    # predict using deepCompare
    pred_ref=compute_predictions(df_copy["seq_ref"].tolist(),device=device)
    pred_alt=compute_predictions(df_copy["seq_alt"].tolist(),device=device)
    
    return pred_alt-pred_ref
    