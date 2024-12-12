import pandas as pd
from loguru import logger
import argparse


from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_index, calculate_tf_cooperativity_index


# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------

def read_all_files(tpm_thresh,nonlinearity_thresh):
    # read cooperativity
    df_promoter_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_hepg2.csv",tpm_thresh,nonlinearity_thresh,track_nums=[0,2,4,6])
    df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_k562.csv",tpm_thresh,nonlinearity_thresh,track_nums=[1,3,5,7])
    df_enhancer_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_hepg2.csv",tpm_thresh,nonlinearity_thresh,track_nums=[0,2,4,6])
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_k562.csv",tpm_thresh,nonlinearity_thresh,track_nums=[1,3,5,7])
    df_promoter_hepg2["dataset"]="promoter_hepg2"
    df_promoter_k562["dataset"]="promoter_k562"
    df_enhancer_hepg2["dataset"]="enhancer_hepg2"
    df_enhancer_k562["dataset"]="enhancer_k562"
    df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562],axis=0)
    return df


def write_ci_by_cell_type(df_orig,cell_type,tpm_thresh,nonlinearity_thresh):
    df=df_orig.copy()
    df=df[df["cell_line"]==cell_type].reset_index(drop=True)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_type}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{tpm_thresh}_{nonlinearity_thresh}_{cell_type}.csv",index=False)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tpm_thresh", type=float)
    parser.add_argument("--nonlinearity_thresh", type=float)
    args = parser.parse_args()

    df=read_all_files(args.tpm_thresh,args.nonlinearity_thresh)
    
    df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[-1])
    for suffix in ["hepg2","k562"]:
        logger.info(f"Processing {suffix}")
        write_ci_by_cell_type(df,suffix,args.tpm_thresh,args.nonlinearity_thresh)



