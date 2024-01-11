import pandas as pd
import sys
import pybedtools
import seaborn as sns
import numpy as np

sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from gradxinp import compute_gradxinp_from_seq
from utils import SeqExtractor
from motif_annotation import JasparAnnotator, ReMapAnnotator, add_feat_imp 

#-----------------------
# Functions
#-----------------------

def add_tf_expr(motif_df,tf_expr):
    """Add TF expression to the motif data frame.
    Args:
        motif_df: A data frame of start, end, TF, score and strand
        tf_expr: A data frame of TF expression with only two columns:'HGNC.symbol' and 'TPM'
    Returns:
        A data frame of start, end, TF, score, strand and TPM
    """
    motif_df=pd.merge(motif_df,tf_expr,left_on='protein',right_on='HGNC.symbol',how='left')
    motif_df=motif_df.loc[:,['chromosome','start','end','strand','protein','score','TPM']]
    # deal with dimer, i.e. X::Y
    # step 1: find all dimers, the "protein" with two colons
    dimers_df=motif_df[motif_df['protein'].str.contains('::')].copy()
    dimers_df.drop(columns='TPM',inplace=True)
    # step 2: split the dimers into two columns
    dimers_df[['tf1','tf2']]=dimers_df['protein'].str.split('::',expand=True)
    # step 3: merge the dimers_df with tf_expr twice, one for tf1 and one for tf2
    dimers_df=pd.merge(dimers_df,tf_expr,left_on='tf1',right_on='HGNC.symbol',how='left')
    dimers_df.drop(columns='HGNC.symbol', inplace=True)
    dimers_df.rename(columns={'TPM':'TPM1'},inplace=True)
    dimers_df=pd.merge(dimers_df,tf_expr,left_on='tf2',right_on='HGNC.symbol',how='left')
    dimers_df.drop(columns='HGNC.symbol', inplace=True)
    dimers_df.rename(columns={'TPM':'TPM2'},inplace=True)
    # step 4: calculate the minimum TPM of the two TFs
    dimers_df['TPM']=dimers_df[['TPM1','TPM2']].min(axis=1)
    # step 5: drop the redundant columns
    dimers_df.drop(columns=['tf1','tf2','TPM1','TPM2'],inplace=True)
    # step 6: merge the dimers_df with motif_df
    motif_df=motif_df[~motif_df['protein'].str.contains('::')].copy()
    motif_df=pd.concat([motif_df,dimers_df],axis=0,ignore_index=True)
    motif_df.sort_values(by=['start','end'],inplace=True)
    return motif_df

def merge_intervals(df, other_cols=['protein', 'max_gradxinp','mean_gradxinp','mean_abs_gradxinp'], operations=["mode","mean","mean","mean"]):
    bed = pybedtools.BedTool.from_dataframe(df)
    col_idxs = [df.columns.get_loc(col) + 1 for col in other_cols]  # +1 because BedTool columns are 1-indexed
    col_str = ','.join(map(str, col_idxs))
    op_str = ','.join(operations)
    merged = bed.merge(c=col_str, o=op_str).to_dataframe(names=['chrom', 'start', 'end'] + other_cols)
    return merged
#-----------------------
# Analysis
#-----------------------
# load data and tools
seq_extractor = SeqExtractor()
jaspar_annotator=JasparAnnotator()
remap_annotator = ReMapAnnotator("/isdata/alab/people/pcr980/Resource/ReMap2022/ReMap2022_K-562.bed")

# Define genomic regions
region = ("chr10", 132396935, 132397534)

# annotate with jaspar
motif_df=jaspar_annotator.annotate(region)


# get sequence from regions
seq_extractor=SeqExtractor()
seq=seq_extractor.get_seq(*region)

# annotate with gradxinp
gradxinp_value=compute_gradxinp_from_seq(seq,targets=1)
motif_df['motif_sequence'] = motif_df.apply(lambda row: seq_extractor.get_seq(row["chromosome"], row["start"], row["end"]), axis=1)
motif_df=add_feat_imp(motif_df,region,gradxinp_value)

# annotate with ReMap
motif_df = remap_annotator.annotate(motif_df,region)



# determine the order of motif mutation.
motif_df_sub=motif_df[motif_df['chip_evidence']==True].copy().reset_index(drop=True)
motif_df_sub=motif_df_sub.groupby('protein').apply(merge_intervals).reset_index(drop=True)

# sort motif_df_sub by max_gradxinp from high to low
motif_df_sub.sort_values(by=['max_gradxinp','mean_abs_gradxinp'],ascending=False,inplace=True)
motif_df_sub.reset_index(drop=True,inplace=True)












#-----------------------
# Archived
#-----------------------
# # add tf expression
# tf_expr=pd.read_csv("/binf-isilon/alab/people/pcr980/DeepCompare/Curate_motif_annot/TF_expression_k562.csv",index_col=0)
# tf_expr=tf_expr.loc[:,['HGNC.symbol','TPM']] # column "HGNC.symbol" contains only capital letters
# motif_df=add_tf_expr(motif_df,tf_expr)
# np.sum(motif_df['TPM']==0) # 363
