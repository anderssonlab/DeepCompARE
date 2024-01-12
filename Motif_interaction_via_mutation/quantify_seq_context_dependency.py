import pandas as pd
import numpy as np
import sys
from pybedtools import BedTool
from loguru import logger

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor,generate_random_seq
from prediction import compute_predictions
from motif_annotation import JasparAnnotator, ReMapAnnotator, subset_df_by_region, add_feat_imp 
from gradxinp import compute_gradxinp_from_seq
#-----------------------
# Functions
#-----------------------

def merge_intervals(df, 
                    other_cols=['protein'],    #['protein', 'max_gradxinp','mean_gradxinp','mean_abs_gradxinp'],
                    operations=["mode"]  #["mode","mean","mean","mean"]
                    ):
    bed = BedTool.from_dataframe(df)
    col_idxs = [df.columns.get_loc(col) + 1 for col in other_cols]  # +1 because BedTool columns are 1-indexed
    col_str = ','.join(map(str, col_idxs))
    op_str = ','.join(operations)
    merged = bed.merge(c=col_str, o=op_str).to_dataframe(names=['chromosome', 'start', 'end'] + other_cols)
    return merged



def scramble_seq_spare_motif(seq, motif_start, motif_end,method):
    """
    Scramble the sequence between multiple motif start and end positions.
    Args:
        seq: A string of sequence.
        motif_start
        motif_end
    Returns:
        A string of scrambled sequence with only motif_start to motif_end remain the same
    """
    if method=="N":
        seq_mut="N"*len(seq)
    elif method=="R":
        seq_mut=generate_random_seq(len(seq))
    else:
        raise ValueError("method must be N or R")
    seq_mut=seq_mut[0:motif_start]+seq[motif_start:(motif_end+1)]+seq_mut[(motif_end+1):]
    return seq_mut

    

def extract_motif_info(region):
    """
    Given a genomic region, extract the motif information.
    Args:
        region: A tuple of (chromosome, start, end).
    Returns:
        A data frame of location, TF, Jaspar score and strand.
    """
    motif_df=jaspar_annotator.annotate(region)
    motif_df = remap_annotator.annotate(motif_df,region)
    motif_df=motif_df[motif_df['chip_evidence']==True].copy().reset_index(drop=True)
    motif_df=motif_df.groupby('protein').apply(merge_intervals).reset_index(drop=True)
    motif_df["start_rel"]=motif_df["start"]-region[1]
    motif_df["end_rel"]=motif_df["end"]-region[1]
    return motif_df


def write_mutation_res_for_one_region(region,out_path,region_idx,method="N"):
    seq = seq_extractor.get_seq(*region)
    df_mutation=extract_motif_info(region)
    if df_mutation.shape[0]==0:
        return
    df_mutation["region_idx"]=f"Region{region_idx}"
    df_mutation["seq_orig"]=seq
    df_mutation['seq_mut'] = df_mutation.apply(lambda row: scramble_seq_spare_motif(seq,row['start_rel'],row['end_rel'], method=method), axis=1)
    df_mutation["pred_orig"]=compute_predictions(df_mutation["seq_orig"])[:,1]
    df_mutation["pred_mut"]=compute_predictions(df_mutation["seq_mut"])[:,1]
    df_feat_imp_orig=compute_gradxinp_from_seq(df_mutation["seq_orig"],targets=1)
    df_feat_imp_mut=compute_gradxinp_from_seq(df_mutation["seq_mut"],targets=1)
    df_mutation["feat_imp_orig"]=[np.max(df_feat_imp_orig.iloc[i, row['start_rel']:row['end_rel']+1]) for i, row in df_mutation.iterrows()]
    df_mutation["feat_imp_mut"]=[np.max(df_feat_imp_mut.iloc[i, row['start_rel']:row['end_rel']+1]) for i, row in df_mutation.iterrows()]
    df_mutation.drop(columns=["seq_orig","seq_mut"],inplace=True)
    if region_idx==0:
        df_mutation.to_csv(out_path,mode="w",header=True,index=False)
        return
    df_mutation.to_csv(out_path,mode="a",header=False,index=False)
    return

#-----------------------
# Analysis
#-----------------------
# load data and tools
seq_extractor = SeqExtractor()
jaspar_annotator=JasparAnnotator()
remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")


df_regions = pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/resize_600bp_CAGE_K562.bed",sep='\t',header=None)
df_regions.iloc[:,2]=df_regions.iloc[:,2]-1

for idx in range(df_regions.shape[0]):
    if idx%50==0:
        logger.info(f"{idx} regions processed.")
    region=df_regions.iloc[idx,[0,1,2]]
    write_mutation_res_for_one_region(region,
                                      f"/isdata/alab/people/pcr980/DeepCompare/Motif_interaction_via_mutation/df_spare_motif_cage_k562.csv",
                                      idx,
                                      method="N")


