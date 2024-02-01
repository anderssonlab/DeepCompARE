import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

# step 1: calculate GC content for each region in /isdata/alab/people/pcr980/DeepCompare/Pd2_sequences/seqs_CAGE_K562.csv
df=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd2_sequences/seqs_CAGE_K562.csv')
df.index = df.index.map(lambda x: 'Seq'+str(x))
# to do: calculate 10, 50, 100, 300 GC content

df["GC_content"] = df["sequence"].apply(lambda x: (x.count("G")+x.count("C"))/len(x))


# step 2: read in /isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/motif_df.csv
motif_df=pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/motif_df.csv') # (63331556, 17)
# merge with GC content
motif_df = motif_df.merge(df[["GC_content"]], left_on="seq_idx", right_index=True) 


# step 3: match for motif score
# to do: find better way to adjust for motif score


# subset motif_df to only include TFs with chip_evidence
motif_df=motif_df[motif_df["chip_evidence"]==True] # (4843557, 18)

df_res=pd.DataFrame(columns=["protein","corr","pval"])

tfs=motif_df["protein"].unique()
for tf in tfs:
    df_sub=motif_df[motif_df["protein"]==tf]
    corr, pval = pearsonr(df_sub["max_gradxinp"], df_sub["GC_content"])
    df_res=df_res.append({"protein":tf,"corr":corr,"pval":pval},ignore_index=True)
    

df_res.to_csv('/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/gc_context_impact_on_featimp.csv',index=False)
