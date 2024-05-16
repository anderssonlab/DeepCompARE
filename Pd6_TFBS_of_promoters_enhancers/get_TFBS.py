#!/usr/bin/env python3

import os
import pyranges as pr
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
from motif_annotation import JasparAnnotator, ReMapAnnotator


#-----------------------
# Functions
#-----------------------


def annotate_one_region(jaspar_annotator,
                        remap_annotator,
                        region,
                        out_path,
                        idx):
    
    motif_df=jaspar_annotator.annotate(region)
    motif_df=remap_annotator.annotate(motif_df,region)
    motif_df=motif_df[motif_df["chip_evidence"]==True].reset_index(drop=True)
    if motif_df.shape[0]==0:
        return
    
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    
    # add sequence information
    motif_df["seq_idx"]= f"Seq{idx}"
    # write to csv
    if not os.path.exists(out_path):
        motif_df.to_csv(out_path,index=False)
    else:
         motif_df.to_csv(out_path,index=False,mode="a",header=False)



#-----------------------
# Analysis
#-----------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate motif info')
    parser.add_argument('--file_name', type=str)
    args=parser.parse_args()
    
    # load promoter/enhancer
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{args.file_name}.bed").df
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    # 
    # load tools
    seq_extractor = SeqExtractor()
    jaspar_annotator=JasparAnnotator()
    if "k562" in args.file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")
    if "hepg2" in args.file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")
        
    
    for idx in range(df_regions.shape[0]):
        if idx%1000==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,0:3].tolist()
        annotate_one_region(region,
                            f"motif_info_{args.file_name}.csv",
                            idx,
                            seq_extractor)
    logger.info(f"Finished annotating {args.file_name}.")

