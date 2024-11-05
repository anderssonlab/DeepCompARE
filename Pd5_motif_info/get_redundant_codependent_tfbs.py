import pandas as pd




def add_tf_codependency(df,suffix):
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_codependent_{suffix}.txt", header=None).iloc[:,0].tolist()
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tfs_redundant_{suffix}.txt", header=None).iloc[:,0].tolist()
    df["cooperativity"]="unknown"
    df.loc[df["protein"].isin(tfs_codependent), "cooperativity"]="codependent"
    df.loc[df["protein"].isin(tfs_redundant), "cooperativity"]="redundant"
    return df



for cell_type in ["hepg2","k562"]:
    df1=pd.read_csv(f'motif_info_thresh_500_promoters_{cell_type}.csv')
    df2=pd.read_csv(f'motif_info_thresh_500_enhancers_{cell_type}.csv')
    df=pd.concat([df1,df2],axis=0).reset_index(drop=True)
    df=add_tf_codependency(df,cell_type)
    # select only cooperativity="codependent" or "redundant"
    df_codependent=df[df["cooperativity"]=="codependent"].reset_index(drop=True)
    df_redundant=df[df["cooperativity"]=="redundant"].reset_index(drop=True)
    # save first 3 columns
    df_codependent.iloc[:,0:3].to_csv(f"tfbs_{cell_type}_codependent.csv",index=False)
    df_redundant.iloc[:,0:3].to_csv(f"tfbs_{cell_type}_redundant.csv",index=False)