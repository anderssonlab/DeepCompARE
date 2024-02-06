import pandas as pd
import sys


sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from utils import SeqExtractor,generate_random_seq
from prediction import compute_predictions
from motif_annotation import JasparAnnotator, ReMapAnnotator, subset_df_by_region, add_feat_imp 
from region_ops import merge_intervals
from loguru import logger
#-----------------------
# Functions
#-----------------------


def find_non_overlapping_combinations(df):
    # Assuming subset_df_by_region is a function that cannot be vectorized
    # Use list comprehension for efficiency
    combination_indices = [
        (idx1, idx2)
        for idx1, row in df.iterrows()
        for idx2 in subset_df_by_region(df, row, by="reverse").index
    ]
    # Convert to a set of frozensets for efficient duplicate removal
    combination_indices = set(frozenset(x) for x in combination_indices)
    # Convert back to list of tuples
    combination_indices = [tuple(x) for x in combination_indices]
    return combination_indices



def test_find_non_overlapping_combinations():
    df=pd.DataFrame([["chr1",1,10,"a"],["chr1",5,15,"b"],["chr1",20,30,"c"]],columns=['chromosome','start','end','protein'])
    res=find_non_overlapping_combinations(df)
    assert set(res)==set([(0,2),(1,2)])


def scramble_motifs(seq, motif_starts, motif_ends,method):
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
        if method=="random":
            motif_scrambled = generate_random_seq(len(motif))
        if method=="N":
            motif_scrambled = "N" * len(motif)  
        # Append non-motif and scrambled motif parts
        seq_scrambled += seq[previous_end:start] + motif_scrambled
        previous_end = end
    # Append the remaining part of the sequence if any
    seq_scrambled += seq[previous_end:]
    return seq_scrambled


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
    motif_df=extract_motif_info(region)
    if motif_df.shape[0]==0:
        return
    combination_indices=find_non_overlapping_combinations(motif_df)
    df_mutation = pd.DataFrame(combination_indices, columns=['idx1', 'idx2'])
    if df_mutation.shape[0]==0:
        return
    df_mutation["region_idx"]=f"Region{region_idx}"
    df_mutation = df_mutation.merge(motif_df, left_on='idx1', right_index=True, suffixes=('', '1'))
    df_mutation = df_mutation.merge(motif_df, left_on='idx2', right_index=True, suffixes=('', '2'))
    df_mutation.rename(columns={'start_rel':'start_rel1','end_rel':'end_rel1','protein':'protein1'},inplace=True) 
    df_mutation["seq_orig"]=seq
    df_mutation['seq_mut1'] = df_mutation.apply(lambda row: scramble_motifs(seq, [row['start_rel1']], [row['end_rel1']], method=method), axis=1)
    df_mutation['seq_mut2'] = df_mutation.apply(lambda row: scramble_motifs(seq, [row['start_rel2']], [row['end_rel2']], method=method), axis=1)
    df_mutation['seq_mut_both'] = df_mutation.apply(lambda row: scramble_motifs(seq, 
                                                                                  [row['start_rel1'], row['start_rel2']], 
                                                                                 [row['end_rel1'], row['end_rel2']], method=method), axis=1)
    df_mutation["pred_orig"]=compute_predictions(df_mutation["seq_orig"])[:,1]
    df_mutation["pred_mut1"]=compute_predictions(df_mutation["seq_mut1"])[:,1]
    df_mutation["pred_mut2"]=compute_predictions(df_mutation["seq_mut2"])[:,1]
    df_mutation["pred_mut_both"]=compute_predictions(df_mutation["seq_mut_both"])[:,1]
    # positive ISM score indicate the original motif is contributing to RE activity.
    df_mutation["ism_score_mut1"]=df_mutation["pred_orig"]-df_mutation["pred_mut1"]
    df_mutation["ism_score_mut2"]=df_mutation["pred_orig"]-df_mutation["pred_mut2"]
    df_mutation["ism_score_mut_both"]=df_mutation["pred_orig"]-df_mutation["pred_mut_both"]
    # perfect AND relationship gives 0. The more negative, the less they depend on each other, or the dependency is biased
    df_mutation["AND_relation"]= (df_mutation["ism_score_mut1"]+df_mutation["ism_score_mut2"])/2 - df_mutation["ism_score_mut_both"] 
    # perfect additivity (independence) gives 0. The more negative, the more synergy
    df_mutation["Additivity"]= df_mutation["ism_score_mut1"]+df_mutation["ism_score_mut2"] - df_mutation["ism_score_mut_both"] 
    # remove the redundant columns
    df_mutation.drop(columns=['seq_mut1','seq_mut2','seq_mut_both','seq_orig'],inplace=True)
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
                                      f"/isdata/alab/people/pcr980/DeepCompare/Motif_interaction_via_mutation/df_mutate_motif_cage_k562.csv",
                                      idx,
                                      method="N")







#-----------------------
# Archived
#-----------------------

# regions = [
#            ("chr10", 132396935, 132397534),
#            ("chr10", 8054655, 8055254),
#            ("chr15", 37098712, 37099311),
#            ("chr12", 132828699, 132829298),
#            ("chr13", 114234732, 114235331),
#            ("chr17", 36544867, 36545466)
#            ]



# # add tf expression

# def add_tf_expr(motif_df,tf_expr):
#     """Add TF expression to the motif data frame.
#     Args:
#         motif_df: A data frame of start, end, TF, score and strand
#         tf_expr: A data frame of TF expression with only two columns:'HGNC.symbol' and 'TPM'
#     Returns:
#         A data frame of start, end, TF, score, strand and TPM
#     """
#     motif_df=pd.merge(motif_df,tf_expr,left_on='protein',right_on='HGNC.symbol',how='left')
#     motif_df=motif_df.loc[:,['chromosome','start','end','strand','protein','score','TPM']]
#     # deal with dimer, i.e. X::Y
#     # step 1: find all dimers, the "protein" with two colons
#     dimers_df=motif_df[motif_df['protein'].str.contains('::')].copy()
#     dimers_df.drop(columns='TPM',inplace=True)
#     # step 2: split the dimers into two columns
#     dimers_df[['tf1','tf2']]=dimers_df['protein'].str.split('::',expand=True)
#     # step 3: merge the dimers_df with tf_expr twice, one for tf1 and one for tf2
#     dimers_df=pd.merge(dimers_df,tf_expr,left_on='tf1',right_on='HGNC.symbol',how='left')
#     dimers_df.drop(columns='HGNC.symbol', inplace=True)
#     dimers_df.rename(columns={'TPM':'TPM1'},inplace=True)
#     dimers_df=pd.merge(dimers_df,tf_expr,left_on='tf2',right_on='HGNC.symbol',how='left')
#     dimers_df.drop(columns='HGNC.symbol', inplace=True)
#     dimers_df.rename(columns={'TPM':'TPM2'},inplace=True)
#     # step 4: calculate the minimum TPM of the two TFs
#     dimers_df['TPM']=dimers_df[['TPM1','TPM2']].min(axis=1)
#     # step 5: drop the redundant columns
#     dimers_df.drop(columns=['tf1','tf2','TPM1','TPM2'],inplace=True)
#     # step 6: merge the dimers_df with motif_df
#     motif_df=motif_df[~motif_df['protein'].str.contains('::')].copy()
#     motif_df=pd.concat([motif_df,dimers_df],axis=0,ignore_index=True)
#     motif_df.sort_values(by=['start','end'],inplace=True)
#     return motif_df

# tf_expr=pd.read_csv("/binf-isilon/alab/people/pcr980/DeepCompare/Curate_motif_annot/TF_expression_k562.csv",index_col=0)
# tf_expr=tf_expr.loc[:,['HGNC.symbol','TPM']] # column "HGNC.symbol" contains only capital letters
# motif_df=add_tf_expr(motif_df,tf_expr)
# np.sum(motif_df['TPM']==0) # 363






# # annotate with gradxinp
# gradxinp_value=compute_gradxinp_from_seq(seq,targets=1)
# motif_df['motif_sequence'] = motif_df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
# motif_df=add_feat_imp(motif_df,region,gradxinp_value)
