import pandas as pd
import argparse
import sys
import os


sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from seq_ops import SeqExtractor
from prediction import compute_predictions
from motif_annotation import JasparAnnotator, ReMapAnnotator
from loguru import logger
#-----------------------
# Functions
#-----------------------


def get_combinations(df):
    n=df.shape[0]
    combinations=[(i,j) for i in range(n) for j in range(i+1,n)]
    return combinations


def scramble_motifs(seq, motif_starts, motif_ends):
    """
    Scramble the sequence between multiple motif start and end positions.
    Args:
        seq: A string of sequence.
        motif_starts: A list of integers for motif starts.
        motif_ends: A list of integers for motif ends.
    Returns:
        A string of scrambled sequence.
    """
    if isinstance(motif_starts, int):
        motif_starts = [motif_starts]
    if isinstance(motif_ends, int):
        motif_ends = [motif_ends]
    if len(motif_starts) != len(motif_ends):
        raise ValueError("motif_starts and motif_ends must have the same length")
    # Sort the motifs by start position
    motifs = sorted(zip(motif_starts, motif_ends), key=lambda x: x[0])
    # Initialize variables
    seq_scrambled = ''
    previous_end = 0
    # Iterate and scramble each motif
    for start, end in motifs:
        if start < previous_end:
            raise ValueError("Overlapping motifs detected")
        end = end + 1  # Adjust the end position
        motif = seq[start:end]
        motif_scrambled = "N" * len(motif)  
        # Append non-motif and scrambled motif parts
        seq_scrambled += seq[previous_end:start] + motif_scrambled
        previous_end = end
    # Append the remaining part of the sequence if any
    seq_scrambled += seq[previous_end:]
    return seq_scrambled


def extract_motif_info(jaspar_annotator,remap_annotator,region):
    """
    Given a genomic region, extract the motif information.
    Args:
        region: A tuple of (chromosome, start, end).
    Returns:
        A data frame of location, TF, Jaspar score and strand.
    """
    motif_df= jaspar_annotator.annotate(region)
    motif_df=motif_df[motif_df["score"]>500].reset_index(drop=True)
    if motif_df.shape[0]==0:
        return motif_df
    motif_df=remap_annotator.annotate(motif_df,region)
    motif_df=motif_df[motif_df['chip_evidence']==True].reset_index(drop=True)
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    motif_df["motif_length"]=motif_df["end_rel"]-motif_df["start_rel"]
    motif_df.drop(columns=['chip_evidence'],inplace=True)
    return motif_df


def write_mutation_res_for_one_region(seq_extractor,
                                      jaspar_annotator,
                                      remap_annotator,
                                      region,
                                      out_path,
                                      region_idx,
                                      track_num,
                                      device):
    seq = seq_extractor.get_seq(*region)
    motif_df=extract_motif_info(jaspar_annotator,remap_annotator,region)
    if motif_df.shape[0]==0:
        return
    combination_indices=get_combinations(motif_df)
    df_mutation = pd.DataFrame(combination_indices, columns=['idx1', 'idx2'])
    df_mutation["region_idx"]=f"Region{region_idx}"
    df_mutation = df_mutation.merge(motif_df, left_on='idx1', right_index=True, suffixes=('', '1'))
    df_mutation = df_mutation.merge(motif_df, left_on='idx2', right_index=True, suffixes=('', '2'))
    df_mutation.rename(columns={'start_rel':'start_rel1','end_rel':'end_rel1','protein':'protein1'},inplace=True) 
    
    # remove overlapping motif locations
    df_mutation = df_mutation[df_mutation['start_rel2']>df_mutation['end_rel1']].reset_index(drop=True)
    if df_mutation.shape[0]==0:
        return
    # get sequences
    df_mutation["seq_orig"]=seq
    df_mutation['seq_mut1'] = df_mutation.apply(lambda row: scramble_motifs(seq, row['start_rel1'], row['end_rel1']), axis=1)
    df_mutation['seq_mut2'] = df_mutation.apply(lambda row: scramble_motifs(seq, row['start_rel2'], row['end_rel2']), axis=1)
    df_mutation['seq_mut_both'] = df_mutation.apply(lambda row: scramble_motifs(seq, 
                                                                                [row['start_rel1'], row['start_rel2']],
                                                                                [row['end_rel1'], row['end_rel2']]), 
                                                    axis=1)
    # get predictions
    df_mutation["pred_orig"]=compute_predictions(df_mutation["seq_orig"].tolist(),device=device)[:,track_num]
    df_mutation["pred_mut1"]=compute_predictions(df_mutation["seq_mut1"].tolist(),device=device)[:,track_num]
    df_mutation["pred_mut2"]=compute_predictions(df_mutation["seq_mut2"].tolist(),device=device)[:,track_num]
    df_mutation["pred_mut_both"]=compute_predictions(df_mutation["seq_mut_both"])[:,track_num]
    # positive ISM score indicate the original motif is contributing positively to RE activity.
    df_mutation["ism_score_mut1"]=df_mutation["pred_orig"]-df_mutation["pred_mut1"]
    df_mutation["ism_score_mut2"]=df_mutation["pred_orig"]-df_mutation["pred_mut2"]
    df_mutation["ism_score_mut_both"]=df_mutation["pred_orig"]-df_mutation["pred_mut_both"]
    
    df_mutation.drop(columns=['seq_mut1','seq_mut2','seq_mut_both','seq_orig'],inplace=True)
    # if out_path does not exist, create it using mode="w", include header
    if not os.path.exists(out_path):
        df_mutation.to_csv(out_path,mode="w",header=True,index=False)
        return
    df_mutation.to_csv(out_path,mode="a",header=False,index=False)
    return 

def full_analysis(file_prefix,remap_cell_type,track_num,device):
    # load data and tools
    seq_extractor = SeqExtractor()
    jaspar_annotator=JasparAnnotator()
    remap_annotator = ReMapAnnotator(f"/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_{remap_cell_type}.bed")
    
    # read the regions
    df_regions = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_prefix}.bed",sep=" ",header=None)

    for idx in range(df_regions.shape[0]):
        if idx%1000==0:
            logger.info(f"{idx} regions processed.")
        region=df_regions.iloc[idx,[0,1,2]]
        write_mutation_res_for_one_region(seq_extractor,
                                          jaspar_annotator,
                                          remap_annotator,
                                          region,
                                          f"df_mutate_pair_{file_prefix}_remap_{remap_cell_type},track{track_num}.csv",
                                          idx,
                                          track_num,
                                          device)



if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--file_prefix",type=str)
    parser.add_argument("--remap_cell_type",type=str)
    parser.add_argument("--track_num",type=int)
    parser.add_argument("--device",type=str)
    args=parser.parse_args()
    
    full_analysis(args.file_prefix,args.remap_cell_type,args.track_num,args.device)
    logger.info(f"Done with {args.file_prefix} {args.remap_cell_type} {args.track_num}")