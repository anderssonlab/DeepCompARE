import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_ratio, get_ppi



# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def read_all_files():
    # read cooperativity
    #df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_promoters_k562.csv")
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/mutate_pairs_enhancers_k562.csv")
    #df=pd.concat([df_promoter_k562,df_enhancer_k562],axis=0)
    df=df_enhancer_k562
    return df




# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------
df=read_all_files()
df=calculate_tf_pair_cooperativity_ratio(df)
df=df[df["sum_ci"]>1].reset_index(drop=True)
df.to_csv("tf_pair_cooperativity_ratio_enhancer_alone.csv",index=False)


df=pd.read_csv("tf_pair_cooperativity_ratio.csv")
get_ppi(df,"redundancy",0.1,"redundancy_ppi.csv")
get_ppi(df,"codependency",0.9,"codependency_ppi.csv")


df=pd.read_csv("tf_pair_cooperativity_ratio_enhancer_alone.csv")
get_ppi(df,"redundancy",0.1,"redundancy_ppi_enhancer_alone.csv")
get_ppi(df,"codependency",0.9,"codependency_ppi_enhancer_alone.csv")





