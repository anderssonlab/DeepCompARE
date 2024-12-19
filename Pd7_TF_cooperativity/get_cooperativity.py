import pandas as pd
import sys

from loguru import logger
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_index, calculate_tf_cooperativity_index


# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------

def read_all_files(threshold):
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
    return df


def write_ci_by_cell_type(df_orig,cell_type):
    df=df_orig.copy()
    df=df[df["cell_line"]==cell_type].reset_index(drop=True)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{cell_type}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{cell_type}.csv",index=False)



def write_ci_by_re(df_orig,re_type):
    df=df_orig.copy()
    df=df[df["re"]==re_type].reset_index(drop=True)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{re_type}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{re_type}.csv",index=False)



def write_tfs_codependent_and_redundant(path,suffix):
    df_tf=pd.read_csv(path)
    df_tf=df_tf[df_tf["c_sum"]>5].reset_index(drop=True)
    # redundant_tfs have cooperativity_index<0.3
    redundant_tfs=df_tf[df_tf["cooperativity_index"]<0.3]["protein2"].to_list()
    # sort alphabetically
    redundant_tfs.sort()
    with open(f"tfs_redundant_{suffix}.txt","w") as f:
        f.write("\n".join(redundant_tfs))
    codependent_tfs=df_tf[df_tf["cooperativity_index"]>0.7]["protein2"].to_list()
    codependent_tfs.sort()
    with open(f"tfs_codependent_{suffix}.txt","w") as f:
        f.write("\n".join(codependent_tfs))




df=read_all_files(0.1)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[-1])
for suffix in ["hepg2","k562"]:
    logger.info(f"Processing {suffix}")
    write_ci_by_cell_type(df,suffix)




df["re"]=["promoter" if "promoter" in x else "enhancer" for x in df["dataset"]]
for suffix in ["promoter","enhancer"]:
    logger.info(f"Processing {suffix}")
    write_ci_by_re(df,suffix)


write_tfs_codependent_and_redundant("tf_cooperativity_index_hepg2.csv","hepg2")
write_tfs_codependent_and_redundant("tf_cooperativity_index_k562.csv","k562")

# nohup python3 get_cooperativity.py > get_cooperativity.out &