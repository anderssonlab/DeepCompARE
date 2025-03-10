import pandas as pd


import sys
from loguru import logger
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")



# TODO: refactor with assign_cooperativity

def write_tfs_by_class(path,suffix,thresh_redun,thresh_codep):
    df_tf=pd.read_csv(path)
    
    # write linear 
    with open(f"tfs_linear_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_linear))
    tfs_redundant=df_nonlinear[df_nonlinear["cooperativity_index"]<thresh_redun]["protein2"].to_list()
    with open(f"tfs_redundant_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_redundant))
    tfs_codependent=df_nonlinear[df_nonlinear["cooperativity_index"]>thresh_codep]["protein2"].to_list()
    with open(f"tfs_codependent_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_codependent))
        
        
# for suffix in ["hepg2_pe","k562_pe","hepg2_dhs","k562_dhs"]:
#     df_coop=pd.read_csv(f"tf_pair_cooperativity_index_{suffix}.csv")
#     df_tf=calculate_tf_cooperativity_index(df_coop)
#     df_tf.to_csv(f"tf_cooperativity_index_{suffix}.csv",index=False)


write_tfs_by_class("tf_cooperativity_index_hepg2_pe.csv","hepg2_pe",0.3,0.7)
write_tfs_by_class("tf_cooperativity_index_k562_pe.csv","k562_pe",0.3,0.7)

# TODO: get new thresholds
# write_tfs_by_class("tf_cooperativity_index_hepg2_dhs.csv","hepg2_dhs",0.5,0.81)
# write_tfs_by_class("tf_cooperativity_index_k562_dhs.csv","k562_dhs",0.46,0.83)



