#!/usr/bin/env python3

import os
import pyranges as pr
import numpy as np
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
from region_ops import calc_gc_context
from prediction import compute_predictions
from gradxinp import compute_gradxinp_from_seq
from motif_annotation import JasparAnnotator, ReMapAnnotator
#-----------------------
# Functions
#-----------------------

def scramble_seq_spare_motif(seq, motif_start, motif_end):
    """
    Scramble the sequence, but spare motif start and end positions.
    Args:
        seq: A string of sequence.
        motif_start
        motif_end
    Returns:
        A string of scrambled sequence with only motif_start to motif_end remain the same
    """
    seq_mut="N"*len(seq)
    seq_mut=seq_mut[0:motif_start]+seq[motif_start:(motif_end+1)]+seq_mut[(motif_end+1):]
    return seq_mut

def spare_seq_scramble_motif(seq, motif_start, motif_end):
    """
    Scramble the sequence between multiple motif start and end positions.
    Args:
        seq: A string of sequence.
        motif_start
        motif_end
    Returns:
        A string of scrambled sequence with only motif_start to motif_end remain the same
    """
    if motif_end==len(seq)-1:
        seq_mut=seq[0:motif_start]+"N"*(motif_end+1-motif_start)
        # assert len(seq_mut)==len(seq)
        return seq_mut
    seq_mut=seq[0:motif_start]+"N"*(motif_end+1-motif_start)+seq[(motif_end+1):]
    # assert len(seq_mut)==len(seq)
    return seq_mut



def annotate_by_removing_context(motif_df,sequence,track_num):
    """add prediction and feat_imp before and after removing context"""
    motif_df['seq_mut'] = motif_df.apply(lambda row: scramble_seq_spare_motif(sequence,row['start_rel'],row['end_rel']), axis=1)
    motif_df["pred_orig"]=compute_predictions(sequence)[0,track_num]
    # to interpret this pred_mut, remember to subtract baseline: N*600
    motif_df["pred_remove_context"]=compute_predictions(motif_df["seq_mut"].tolist())[:,track_num]
    df_feat_imp_orig=compute_gradxinp_from_seq(sequence,targets=track_num)
    df_feat_imp_mut=compute_gradxinp_from_seq(motif_df["seq_mut"].tolist(),targets=track_num)
    # for debugging
    # print(df_feat_imp_orig.shape)
    # print(df_feat_imp_mut.shape)
    # print(motif_df.shape)
    # for i,row in motif_df.iterrows():
    #     print(i,row["start_rel"],row["end_rel"])
    #     print(df_feat_imp_mut.iloc[i,row["start_rel"]:row["end_rel"]+1])
    motif_df["feat_imp_orig"]=motif_df.apply(lambda row:np.max(df_feat_imp_orig.iloc[0,row["start_rel"]:row["end_rel"]+1]),axis=1)
    motif_df["feat_imp_remove_context"]=[np.max(df_feat_imp_mut.iloc[i,row["start_rel"]:row["end_rel"]+1]) for i,row in motif_df.iterrows()]
    return motif_df


def calc_ism_motif(motif_df,sequence,track_num):
    """add prediction before and after removing context"""
    motif_df['seq_mut'] = motif_df.apply(lambda row:spare_seq_scramble_motif(sequence,row['start_rel'],row['end_rel']), axis=1)
    motif_df["pred_mut_motif"]=compute_predictions(motif_df["seq_mut"])[:,track_num]
    motif_df["ism_motif"]=motif_df["pred_orig"]-motif_df["pred_mut_motif"]
    motif_df.drop(columns=["seq_mut","pred_mut_motif"],inplace=True)
    return motif_df



def annotate_homotypic_clusters(motif_df,column_suffix):
    protein_counts=motif_df["protein"].value_counts()
    motif_df[f"count_all_TFs_{column_suffix}"]=protein_counts.sum()
    # change protein counts to a dataframe
    protein_counts=protein_counts.reset_index()
    motif_df=motif_df.merge(protein_counts,how="left",left_on="protein",right_on="index")
    motif_df.rename(columns={"protein_y":f"count_TF_{column_suffix}",
                             "protein_x":"protein"},inplace=True)
    motif_df.drop(columns=["index"],inplace=True)
    return motif_df


# TODO: change threshold here
def annotate_one_region(region,out_path,idx,seq_extractor,track_num,score_thresh=500):
    
    # add motif location and TF
    motif_df=jaspar_annotator.annotate(region)
    if motif_df.shape[0]==0:
        return
    
    # add homotypic cluster info to motif_df before jaspar score filtering
    motif_df=annotate_homotypic_clusters(motif_df,"no_thresh")

    
    if score_thresh:
        motif_df=motif_df[motif_df["score"]>score_thresh].reset_index(drop=True)
        if motif_df.shape[0]==0:
            return
    
    logger.info(f"Annotating {motif_df.shape[0]} motifs in region {idx}.")
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    
    # add motif sequence
    motif_df['motif_sequence'] = motif_df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
    
    # add whether it has remap evidence
    motif_df = remap_annotator.annotate(motif_df,region)
    sequence=seq_extractor.get_seq(region[0],region[1],region[2])
    
    # add gc context
    motif_df["seq_gc"]= (sequence.count("G")+sequence.count("C"))/len(sequence)
    motif_df["context_gc_2bp"]=calc_gc_context(motif_df,2,seq_extractor)
    motif_df["context_gc_10bp"]=calc_gc_context(motif_df,10,seq_extractor)
    motif_df["context_gc_50bp"]=calc_gc_context(motif_df,50,seq_extractor)
    motif_df["context_gc_100bp"]=calc_gc_context(motif_df,100,seq_extractor)
    
    # add information before and after removing context
    motif_df=annotate_by_removing_context(motif_df,sequence,track_num)
    motif_df=calc_ism_motif(motif_df,sequence,track_num)
    
    # add whether the region overlaps homotypic clusters
    # TODO: change threshold
    motif_df=annotate_homotypic_clusters(motif_df, "thresh_500")
    
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
    parser.add_argument('--file_name', type=str, help='File name to annotate')
    parser.add_argument('--track_num', type=str, help='File name to annotate')
    args=parser.parse_args()
    
    # Load data and tools
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{args.file_name}.bed").df
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    seq_extractor = SeqExtractor()
    jaspar_annotator=JasparAnnotator()

    if "k562" in args.file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")
    if "hepg2" in args.file_name:
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_Hep-G2.bed")

    for idx in range(df_regions.shape[0]):
    # debug purpose
    # for idx in range(2):
        if idx%100==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,:]
        # TODO: change file name relative to "thresh" here
        annotate_one_region(region,f"motif_info_thresh_500_{args.file_name}.csv",idx,seq_extractor,int(args.track_num))
    logger.info(f"Finished annotating {args.file_name}.")

