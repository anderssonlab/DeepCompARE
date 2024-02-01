#!/usr/bin/env python3

import pandas as pd
import sys
from loguru import logger

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor
from gradxinp import read_gradxinp
from motif_annotation import JasparAnnotator, ReMapAnnotator, add_feat_imp
#-----------------------
# Functions
#-----------------------

def annotate_one_region(region,gradxinp_value,out_path,idx):
    motif_df=jaspar_annotator.annotate(region)
    motif_df['motif_sequence'] = motif_df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
    motif_df=add_feat_imp(motif_df,region,gradxinp_value)
    motif_df = remap_annotator.annotate(motif_df,region)
    motif_df["seq_idx"]= f"Seq{idx}"
    if idx==0:
        motif_df.to_csv(out_path,index=False,mode="a",header=True)
    else:
         motif_df.to_csv(out_path,index=False,mode="a",header=False)



#-----------------------
# Analysis
#-----------------------

# Load data and tools
df_regions = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/resize_600bp_CAGE_K562.bed",sep='\t',header=None)
df_gradxinp = read_gradxinp("/isdata/alab/people/pcr980/DeepCompare/Pd2_metadata_and_featimp/gradxinp_CAGE_K562.csv",track_num=1)
seq_extractor = SeqExtractor()
jaspar_annotator=JasparAnnotator()
remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")

assert df_regions.shape[0]==df_gradxinp.shape[0]
for idx in range(df_regions.shape[0]):
    if idx%100==0:
        logger.info(f"{idx} regions processed.")
    region=df_regions.iloc[idx,:]
    gradxinp_value=df_gradxinp.iloc[idx,:]
    annotate_one_region(region,gradxinp_value,"/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/motif_df.csv",idx)

