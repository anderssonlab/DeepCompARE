import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from region_ops import calc_gc_context
from utils import SeqExtractor


# to do: 

seq_extractor = SeqExtractor()
motif_df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif_highlight_chip/motif_df.csv")
motif_df=motif_df[motif_df["chip_evidence"]==True] 
motif_df["context_gc_0bp"]=calc_gc_context(motif_df,0,seq_extractor)
motif_df["context_gc_2bp"]=calc_gc_context(motif_df,2,seq_extractor)
motif_df["context_gc_10bp"]=calc_gc_context(motif_df,10,seq_extractor)
motif_df["context_gc_50bp"]=calc_gc_context(motif_df,50,seq_extractor)
motif_df["context_gc_100bp"]=calc_gc_context(motif_df,100,seq_extractor)
motif_df["context_gc_300bp"]=calc_gc_context(motif_df,300,seq_extractor)



tfs=motif_df["protein"].unique()
df_res=pd.DataFrame(columns=["protein","corr","pval","context_width","feat_imp"])

for tf in tfs:
    df_sub=motif_df[motif_df["protein"]==tf]
    for context_width in [0,2,10,50,100,300]:
        for feat_imp in ["max_gradxinp","mean_abs_gradxinp"]:
            corr,pval=pearsonr(df_sub[f"context_gc_{context_width}bp"],df_sub[feat_imp])
            df_res=df_res.append({"protein":tf,"corr":corr,"pval":pval,"context_width":context_width,"feat_imp":feat_imp},ignore_index=True)

df_res.to_csv('/isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/corr_gc_context_and_featimp.csv',index=False)
