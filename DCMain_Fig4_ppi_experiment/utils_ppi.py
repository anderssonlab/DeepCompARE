import pandas as pd
from scipy.stats import mannwhitneyu



def get_htfs_list():
    # read human transcription factors
    df_htfs=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv",index_col=0)
    # subset for "Is TF?"=="Yes" or "HGNC symbol"=="SMAD2"
    df_htfs=df_htfs[(df_htfs["Is TF?"]=="Yes")|(df_htfs["HGNC symbol"]=="SMAD2")].reset_index(drop=True)
    return df_htfs["HGNC symbol"].tolist()
    




def read_pooled_found_tf_pairs(significant_flag=True):
    # pool
    df_found_pooled=pd.DataFrame()
    for bait in ["BACH1","MAFG","IKZF1","RREB1","RFX5"]:
        df_bait=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
        df_found=df_bait[df_bait["significant"]==significant_flag].reset_index(drop=True)
        df_found["bait"]=bait
        df_found_pooled=pd.concat([df_found_pooled,df_found],axis=0)
    # read human transcription factors
    htfs = get_htfs_list()
    # subset to contain only TFs
    df_found_pooled=df_found_pooled[df_found_pooled["gene"].isin(htfs)].reset_index(drop=True)
    # rename df_found_pooled
    df_found_pooled.rename(columns={"gene":"protein1","bait":"protein2"},inplace=True)
    # reverse protein1 and protein2
    df_found_pooled2=df_found_pooled.copy()
    df_found_pooled2.rename(columns={"protein1":"protein2","protein2":"protein1"},inplace=True)
    df_found_pooled=pd.concat([df_found_pooled,df_found_pooled2],axis=0).reset_index(drop=True)
    # remove protein1==protein2
    df_found_pooled=df_found_pooled[df_found_pooled["protein1"]!=df_found_pooled["protein2"]].reset_index(drop=True)
    # select protein2 in ["BACH1","MAFG","IKZF1","RREB1","RFX5"]
    df_found_pooled=df_found_pooled[df_found_pooled["protein2"].isin(["BACH1","MAFG","IKZF1","RREB1","RFX5"])].reset_index(drop=True)
    # remove duplicates based on protein1 and protein2
    df_found_pooled=df_found_pooled.drop_duplicates(subset=["protein1","protein2"]).reset_index(drop=True)
    df_found_pooled=df_found_pooled[["protein1","protein2"]]
    df_found_pooled["experiment_state"]="significant"
    return df_found_pooled






def read_filtered_file(bait):
    df=pd.read_csv(f"Pd1_PPI_experiment/250107.{bait}.K562.Taplin.FilteredProteins.txt", sep='\t')
    df["reference"]=df["reference"].apply(lambda x: x.split("|")[-1])
    df["gene"]=df["reference"].apply(lambda x: x.split("_")[0])
    df["species"]=df["reference"].apply(lambda x: x.split("_")[-1])
    # remove nonhuman
    df=df[df["species"]=="HUMAN"].reset_index(drop=True)
    df["bait"]=bait
    # retain columns "gene" and "bait", rename to "protein1" and "protein2"
    df=df[["gene","bait"]].rename(columns={"gene":"protein1","bait":"protein2"})
    return df




def read_pooled_subthresh_tf_pairs():
    df_subthresh_pooled=pd.DataFrame()
    for bait in ["BACH1","MAFG","IKZF1","RREB1","RFX5"]:
        df_subthresh=read_filtered_file(bait)
        df_subthresh_pooled=pd.concat([df_subthresh_pooled,df_subthresh],axis=0)
    # read human transcription factors
    htfs = get_htfs_list()
    # subset to contain only TFs
    df_subthresh_pooled=df_subthresh_pooled[df_subthresh_pooled["protein1"].isin(htfs)].reset_index(drop=True)
    # reverse protein1 and protein2
    df_subthresh_pooled2=df_subthresh_pooled.copy()
    df_subthresh_pooled2.rename(columns={"protein1":"protein2","protein2":"protein1"},inplace=True)
    df_subthresh_pooled=pd.concat([df_subthresh_pooled,df_subthresh_pooled2],axis=0).reset_index(drop=True)
    # merge with found but insignificant TFs
    df_found_pooled=read_pooled_found_tf_pairs(significant_flag=False)
    # retain columns protein1, protein2
    df_found_pooled=df_found_pooled[["protein1","protein2"]]
    df_subthresh_pooled=pd.concat([df_subthresh_pooled,df_found_pooled],axis=0).reset_index(drop=True)
    # subset for protein2 in ["BACH1","MAFG","IKZF1","RREB1","RFX5"]
    df_subthresh_pooled=df_subthresh_pooled[df_subthresh_pooled["protein2"].isin(["BACH1","MAFG","IKZF1","RREB1","RFX5"])].reset_index(drop=True)
    # remove duplicates based on protein1 and protein2
    df_subthresh_pooled=df_subthresh_pooled.drop_duplicates(subset=["protein1","protein2"]).reset_index(drop=True)
    df_subthresh_pooled["experiment_state"]="subthreshold"
    return df_subthresh_pooled

        










def mannwhitneyu_with_nan(list1,list2):
    # remove nan
    list1=list1[~pd.isnull(list1)]
    list2=list2[~pd.isnull(list2)]
    stat,p=mannwhitneyu(list1,list2)
    return stat,p

