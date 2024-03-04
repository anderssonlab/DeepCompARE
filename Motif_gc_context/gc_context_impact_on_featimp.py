import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python")
from region_ops import calc_gc_context
from seq_ops import SeqExtractor
from stat_tests import bin_and_label




#-----------------------------------------
# Step 1: get mean imp for each gc content
#-----------------------------------------

seq_extractor = SeqExtractor()
motif_df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif_highlight_chip/motif_df.csv")
# to do: Maybe shouldn't filter out the non-chip evidence motifs
motif_df=motif_df[motif_df["chip_evidence"]==True] 
motif_df["context_gc_0bp"]=calc_gc_context(motif_df,0,seq_extractor)
motif_df["context_gc_2bp"]=calc_gc_context(motif_df,2,seq_extractor)
motif_df["context_gc_10bp"]=calc_gc_context(motif_df,10,seq_extractor)
motif_df["context_gc_50bp"]=calc_gc_context(motif_df,50,seq_extractor)
motif_df["context_gc_100bp"]=calc_gc_context(motif_df,100,seq_extractor)
motif_df["context_gc_300bp"]=calc_gc_context(motif_df,300,seq_extractor)



tfs=motif_df["protein"].unique()


df_res=pd.DataFrame()

for tf in tfs:
    logger.info(f"Processing {tf}")
    df_sub=motif_df[motif_df["protein"]==tf]
    # create bin edges, from 0 to 1 in 0.05 increments
    bins=[i/20 for i in range(21)]
    df_sub=bin_and_label(df_sub,"context_gc_10bp",bin_edges=bins)
    # calculate mean feat_imp for each bin, convert to data frame
    bin_summary=df_sub.groupby("Bin")["max_gradxinp"].median()
    df_bin_summary=pd.DataFrame(bin_summary)
    df_bin_summary["protein"]=tf
    df_res=df_res.append(df_bin_summary)
    

df_res.to_csv('/isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/binned_median_10.csv')


#-----------------------------------------
# Step2: get most preferred gc content for each tf
#-----------------------------------------

# for each tf, find the gc content with the highest max_gradxinp
df_res=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Motif_gc_context/binned_median.csv')
df_gc_preference=pd.DataFrame(df_res.groupby("protein")["max_gradxinp"].idxmax())

df_gc_preference.sort_values(by="max_gradxinp",ascending=True,inplace=True)
