#!/usr/bin/python3
import numpy as np
from prediction import compute_predictions
from region_ops import resize_df
from loguru import logger

def predict_vcf(df_vcf,seq_extractor):
    df_copy=df_vcf.copy()
    logger.info(f"df_copy.shape={df_copy.shape}")
    df_copy=resize_df(df_copy,600)
    df_copy["seq_ref"]=df_copy.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]-1), axis=1)
    logger.info("Reference sequence extracted!")
    # get the center base, validate the center base is same as REF
    df_copy["seq_ref_center"]=df_copy["seq_ref"].apply(lambda x: x[300])
    assert np.all(df_copy["REF"]==df_copy["seq_ref_center"])
    logger.info("Reference bases validated!")
    
    # change center base to ALT
    df_copy["prefix"]=df_copy["seq_ref"].apply(lambda x: x[0:300])
    df_copy["suffix"]=df_copy["seq_ref"].apply(lambda x: x[301:600])
    df_copy["seq_alt"]=df_copy["prefix"]+df_copy["ALT"]+df_copy["suffix"]   
    logger.info("Alternative sequence generated!")
    # validate seq_alt and seq_ref differ only at center base
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[300])!=df_copy["seq_ref"].apply(lambda x: x[300]))
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[299])==df_copy["seq_ref"].apply(lambda x: x[299]))
    assert np.all(df_copy["seq_alt"].apply(lambda x: x[301])==df_copy["seq_ref"].apply(lambda x: x[301]))
    logger.info("Alternative bases validated!")
    
    # predict ref
    pred_ref=compute_predictions(df_copy["seq_ref"].tolist())
    logger.info("Reference sequence predicted!")
    # predict alt
    pred_alt=compute_predictions(df_copy["seq_alt"].tolist())
    logger.info("Alternative sequence predicted!")
    
    return pred_alt-pred_ref
    