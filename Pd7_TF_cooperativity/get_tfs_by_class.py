import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity 



def write_tfs_by_class(path,suffix,thresh_redun,thresh_codep):
    df_tf=pd.read_csv(path)
    df_tf=assign_cooperativity(df_tf,5,0.95,thresh_redun,thresh_codep)
    tfs_linear=df_tf[df_tf["cooperativity"]=="Linear"]["protein2"].tolist()
    tfs_redundant=df_tf[df_tf["cooperativity"]=="Redundant"]["protein2"].tolist()
    tfs_intermediate=df_tf[df_tf["cooperativity"]=="Intermediate"]["protein2"].tolist()
    tfs_codependent=df_tf[df_tf["cooperativity"]=="Codependent"]["protein2"].tolist()
    with open(f"tfs_linear_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_linear))
    with open(f"tfs_redundant_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_redundant))
    with open(f"tfs_intermediate_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_intermediate))
    with open(f"tfs_codependent_{suffix}.txt","w") as f:
        f.write("\n".join(tfs_codependent))


write_tfs_by_class("tf_cooperativity_index_hepg2_pe.csv","hepg2_pe",0.3,0.7)
write_tfs_by_class("tf_cooperativity_index_k562_pe.csv","k562_pe",0.3,0.7)

write_tfs_by_class("tf_cooperativity_index_hepg2_dhs.csv","hepg2_dhs",0.48,0.76)
write_tfs_by_class("tf_cooperativity_index_k562_dhs.csv","k562_dhs",0.44,0.81)



