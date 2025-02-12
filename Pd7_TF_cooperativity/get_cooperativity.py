import pandas as pd
import sys

from loguru import logger
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_index, calculate_tf_cooperativity_index


# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------

def read_all_pe_files(threshold):
    # read cooperativity
    df_promoter_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_enhancer_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_promoter_hepg2["dataset"]="promoter_hepg2"
    df_promoter_k562["dataset"]="promoter_k562"
    df_enhancer_hepg2["dataset"]="enhancer_hepg2"
    df_enhancer_k562["dataset"]="enhancer_k562"
    df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562],axis=0)
    df["re"]=["promoter" if "promoter" in x else "enhancer" for x in df["dataset"]]
    df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[-1])
    return df




def read_all_dhs_files(threshold):
    # read cooperativity
    df_proximal_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_dhs_proximal_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_proximal_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_dhs_proximal_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_distal_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_dhs_distal_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_distal_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_dhs_distal_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_proximal_hepg2["dataset"]="proximal_hepg2"
    df_proximal_k562["dataset"]="proximal_k562"
    df_distal_hepg2["dataset"]="distal_hepg2"
    df_distal_k562["dataset"]="distal_k562"
    df=pd.concat([df_proximal_hepg2,df_proximal_k562,df_distal_hepg2,df_distal_k562],axis=0)
    df["re"]=["proximal" if "proximal" in x else "distal" for x in df["dataset"]]
    df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[-1])
    return df



def write_ci(df,suffix):
    df_coop=calculate_tf_pair_cooperativity_index(df)
    df_coop.to_csv(f"tf_pair_cooperativity_index_{suffix}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df_coop)
    df_tf.to_csv(f"tf_cooperativity_index_{suffix}.csv",index=False)




df=read_all_pe_files(0.1)
for cell in ["hepg2","k562"]:
    logger.info(f"Processing {cell}")
    df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
    write_ci(df_sub,f"{cell}_pe")


for re in ["promoter","enhancer"]:
    logger.info(f"Processing {re} ")
    df_sub=df[df["re"]==re].reset_index(drop=True)
    write_ci(df_sub,f"{cell}_pe")


df=read_all_dhs_files(0.1)
for cell in ["hepg2","k562"]:
    logger.info(f"Processing {cell}")
    df_sub=df[df["cell_line"]==cell].reset_index(drop=True)
    write_ci(df_sub,f"{cell}_dhs")




# nohup python3 get_cooperativity.py > get_cooperativity.out &
