#!/usr/bin/env python3

import os
import pyranges as pr
import numpy as np
import sys
from loguru import logger
import argparse

sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor, shuffle_dinucleotides
from region_ops import calc_gc_context
from prediction import compute_predictions
from gradxinp import compute_gradxinp_from_seq
from seq_annotators import JasparAnnotator, ReMapAnnotator
#-----------------------
# Functions
#-----------------------

def ablate_seq_spare_motif(seq, motif_start, motif_end):
    seq_mut="N"*len(seq)
    return seq_mut[0:motif_start]+seq[motif_start:(motif_end+1)]+seq_mut[(motif_end+1):]

def dishuffle_seq_spare_motif(seq, motif_start, motif_end):
    seq_mut=shuffle_dinucleotides(seq)
    return seq_mut[0:motif_start]+seq[motif_start:(motif_end+1)]+seq_mut[(motif_end+1):]

def spare_seq_ablate_motif(seq, motif_start, motif_end):
    if motif_end==len(seq)-1:
        seq_mut=seq[0:motif_start]+"N"*(motif_end+1-motif_start)
        return seq_mut
    seq_mut=seq[0:motif_start]+"N"*(motif_end+1-motif_start)+seq[(motif_end+1):]
    return seq_mut

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

def extract_motif_gradxinp(motif_df,df_feat_imp):
    return [np.max(df_feat_imp.iloc[0,row["start_rel"]:row["end_rel"]+1]) for _,row in motif_df.iterrows()]

def annotate(motif_df,sequence,track_num,func,column_suffix,device):
    """
    calculate predictions and feat_imp_motif for sequence mutated by func
    """
    motif_df[f"seq_{column_suffix}"] = motif_df.apply(lambda row: func(sequence,row['start_rel'],row['end_rel']), axis=1)
    motif_df[f"pred_{column_suffix}"]=compute_predictions(motif_df[f"seq_{column_suffix}"],device=device)[:,track_num]
    df_feat_imp_mut=compute_gradxinp_from_seq(motif_df[f"seq_{column_suffix}"],targets=track_num,device=device)
    motif_df[f"feat_imp_{column_suffix}"]=extract_motif_gradxinp(motif_df,df_feat_imp_mut)
    motif_df.drop(columns=[f"seq_{column_suffix}"],inplace=True)
    return motif_df


# TODO: change threshold here
def annotate_one_region(region,out_path,idx,seq_extractor,track_num,device,score_thresh=500):
    jaspar_annotator=JasparAnnotator("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2022_hg38.bb")
    motif_df=jaspar_annotator.annotate(region,by="contained")
    if motif_df.shape[0]==0:
        return
    #motif_df=annotate_homotypic_clusters(motif_df,"no_thresh")
    if score_thresh:
        motif_df=motif_df[motif_df["score"]>score_thresh].reset_index(drop=True)
        if motif_df.shape[0]==0:
            return
    
    logger.info(f"Annotating {motif_df.shape[0]} motifs in region {idx}.")
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    
    # add whether it has remap evidence
    motif_df = remap_annotator.annotate(motif_df,region)
    
    # add sequence and motif sequence
    # motif_df['motif_sequence'] = motif_df.apply(lambda row: seq_extractor.get_seq(row[0:3].tolist()), axis=1)
    sequence=seq_extractor.get_seq(region)

    # add gc context
    # motif_df["seq_gc"]= (sequence.count("G")+sequence.count("C"))/len(sequence)
    # motif_df["context_gc_2bp"]=calc_gc_context(motif_df,2,seq_extractor)
    # motif_df["context_gc_10bp"]=calc_gc_context(motif_df,10,seq_extractor)
    # motif_df["context_gc_50bp"]=calc_gc_context(motif_df,50,seq_extractor)
    # motif_df["context_gc_100bp"]=calc_gc_context(motif_df,100,seq_extractor)
    
    # add information before and after removing context
    # motif_df=annotate(motif_df,sequence,track_num,ablate_seq_spare_motif,"ablate_context",device=device)
    motif_df=annotate(motif_df,sequence,track_num,spare_seq_ablate_motif,"ablate_motif",device=device)
    # motif_df=annotate(motif_df,sequence,track_num,dishuffle_seq_spare_motif,"dishuffle_context",device=device)
    # motif_df=annotate_homotypic_clusters(motif_df, "thresh_500")
    
    # add prediction and feat_imp of original sequence
    motif_df["pred_orig"]=compute_predictions(sequence,device=device)[0,track_num]
    # df_feat_imp_orig=compute_gradxinp_from_seq(sequence,targets=track_num,device=device)
    # Note! feat_imp is gradient*input
    # motif_df["feat_imp_orig"]=extract_motif_gradxinp(motif_df,df_feat_imp_orig)
    motif_df["ism_motif"]=motif_df["pred_orig"]-motif_df["pred_ablate_motif"]
    
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
    parser.add_argument('--track_num', type=str)
    parser.add_argument('--device', type=str)
    
    args=parser.parse_args()
    
    # Load data and tools
    #df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{args.file_name}.bed").df
    df_regions=pr.read_bed(f"/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/{args.file_name}.bed").df
    df_regions.iloc[:,2]=df_regions.iloc[:,2]-1
    seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    

    if "k562" in args.file_name.lower():
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_K-562.bed")
    if "hepg2" in args.file_name.lower():
        remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_hg38_Hep-G2.bed")

    for idx in range(df_regions.shape[0]):
        if idx%1000==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,0:3].tolist()
        # TODO: change file name relative to "thresh" here
        annotate_one_region(region,
                            f"motif_info_thresh_500_{args.file_name}_track{args.track_num}.csv",
                            idx,
                            seq_extractor,
                            int(args.track_num),
                            args.device)
    logger.info(f"Finished annotating {args.file_name}, track{args.track_num}.")

