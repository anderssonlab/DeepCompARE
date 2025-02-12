import pandas as pd

def write_tfs_codependent_and_redundant(path,suffix,thresh_redun, thresh_codep):
    df_tf=pd.read_csv(path)
    df_tf=df_tf[df_tf["c_sum"]>1].reset_index(drop=True)
    # redundant_tfs have cooperativity_index<thresh_redun
    redundant_tfs=df_tf[df_tf["cooperativity_index"]<thresh_redun]["protein2"].to_list()
    # sort alphabetically
    redundant_tfs.sort()
    with open(f"tfs_redundant_{suffix}.txt","w") as f:
        f.write("\n".join(redundant_tfs))
    codependent_tfs=df_tf[df_tf["cooperativity_index"]>thresh_codep]["protein2"].to_list()
    codependent_tfs.sort()
    with open(f"tfs_codependent_{suffix}.txt","w") as f:
        f.write("\n".join(codependent_tfs))



write_tfs_codependent_and_redundant("tf_cooperativity_index_hepg2_pe.csv","hepg2_pe",0.3,0.7)
write_tfs_codependent_and_redundant("tf_cooperativity_index_k562_pe.csv","k562_pe",0.3,0.7)


write_tfs_codependent_and_redundant("tf_cooperativity_index_hepg2_dhs.csv","hepg2_dhs",0.5,0.81)
write_tfs_codependent_and_redundant("tf_cooperativity_index_k562_dhs.csv","k562_dhs",0.46,0.83)



