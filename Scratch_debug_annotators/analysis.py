import pandas as pd




def read_file_remove_confusion(file_name,threshold,track_nums=[1,3,5,7]):
    """
    Read in file, remove mistaking rows, calculate cooperativity, and return the dataframe with cooperativity information
    """
    # read file
    df=pd.read_csv(file_name)
    
    # remove rows if protein1 or protein2 are not expressed
    if "hepg2" in file_name:
        expressed_protein_list=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_hepg2.tsv", sep='\t', header=None).iloc[:,0].tolist()
    elif "k562" in file_name:
        expressed_protein_list=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/RNA_expression/expressed_tf_list_k562.tsv", sep='\t', header=None).iloc[:,0].tolist()
    df=df[df["protein1"].isin(expressed_protein_list)].reset_index(drop=True)
    df=df[df["protein2"].isin(expressed_protein_list)].reset_index(drop=True)
    
    tracks_to_remove=[i for i in range(16) if i not in track_nums]
    # remove all columns with track{i} and i is in tracks_to_remove
    col_suffix_to_remove=[f"_track{i}" for i in tracks_to_remove]
    col_names_to_remove=[col for col in df.columns if any([col.endswith(suffix) for suffix in col_suffix_to_remove])]
    df.drop(columns=col_names_to_remove,inplace=True)
    # calculate cooperativity index
    for i in track_nums:
        df[f"isa2_wo_protein1_track{i}"]=df[f'isa_both_track{i}']-df[f'isa1_track{i}']
        # remove potentially mistaking rows and conditional repressors
        df.loc[(df[f"isa1_track{i}"]<0) | (df[f"isa2_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df.loc[(df[f"isa2_wo_protein1_track{i}"]<0),[col for col in df.columns if col.endswith(f"track{i}")]]=None
        df[f"c_track{i}"]=df[f"isa2_track{i}"]-df[f"isa2_wo_protein1_track{i}"]
        # assign cooperativity
        df[f"cooperativity_track{i}"]="unknown"
        df.loc[df[f"c_track{i}"]< -threshold,f"cooperativity_track{i}"]="redundancy"
        df.loc[df[f"c_track{i}"]> threshold,f"cooperativity_track{i}"]="codependency"
    # drop columns starting with "pred" and isa, and motif_length
    col_names_to_remove=[col for col in df.columns if col.startswith("pred") or col.startswith("motif_length")]
    df.drop(columns=col_names_to_remove,inplace=True)
    return df



df_later_subset=read_file_remove_confusion("mutate_pairs_enhancers_k562_head_no_subset.csv",0.1)
df_later_subset=df_later_subset[df_later_subset['region_idx']=="Region0"].reset_index(drop=True)

df_jaspar_subset=read_file_remove_confusion("mutate_pairs_enhancers_k562_head_subset_by_jaspar.csv",0.1)
df_jaspar_subset=df_jaspar_subset[df_jaspar_subset['region_idx']=="Region0"].reset_index(drop=True)


df_later_subset[["protein1","protein2"]]

df_jaspar_subset[["protein1","protein2"]]