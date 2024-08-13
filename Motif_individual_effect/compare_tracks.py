import pandas as pd

df=pd.read_csv("tf_individual_effect_by_cell_type.csv")

def get_cage_vs_sure(df,cell_type):
    df_sub=df[df["dataset"].str.contains(cell_type)].reset_index(drop=True)
    df_sure=df_sub[(df_sub["dstat_ism_cage"]<0) & (df_sub["dstat_ism_sure"] > 0) & (abs(df_sub["dstat_ism_cage"] - df_sub["dstat_ism_sure"]) > 0.3)]
    df_cage=df_sub[(df_sub["dstat_ism_cage"]>0) & (df_sub["dstat_ism_sure"] < 0) & (abs(df_sub["dstat_ism_cage"] - df_sub["dstat_ism_sure"]) > 0.3)]
    tfs_cage=df_cage["protein"].tolist()
    tfs_sure=df_sure["protein"].tolist()
    return tfs_cage,tfs_sure



tfs_cage_hepg2,tfs_sure_hepg2=get_cage_vs_sure(df,"hepg2")
tfs_cage_k562,tfs_sure_k562=get_cage_vs_sure(df,"k562")

# count overlaps
print(len(set(tfs_cage_hepg2) & set(tfs_cage_k562))/len(set(tfs_cage_hepg2) | set(tfs_cage_k562)))
print(len(set(tfs_sure_hepg2) & set(tfs_sure_k562))/len(set(tfs_sure_hepg2) | set(tfs_sure_k562)))


set(tfs_cage_hepg2) & set(tfs_cage_k562)