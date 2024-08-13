import pandas as pd
from scipy.stats import ks_2samp
from loguru import logger

prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_thresh_500_"


#------------------
# Helper functions
#------------------

def read_data(re_type):
    df_list=[]
    if "hepg2" in re_type:
        track_list=[0,2,4,6]
    else:
        track_list=[1,3,5,7]
    for track in track_list:
        df=pd.read_csv(f"{prefix}{re_type}_track{track}.csv")
        df["track"]=f"track{track}"
        df_list.append(df)
    df_res=pd.concat(df_list).reset_index(drop=True)
    df_res=df_res[df_res["chip_evidence"]==True].reset_index(drop=True)
    return df_res





#----------------------------
# Analysis1: write tf individual effect (chip True)
#----------------------------

df_promoters_hepg2=read_data("promoters_hepg2")
df_promoters_k562=read_data("promoters_k562")
df_enhancers_hepg2=read_data("enhancers_hepg2")
df_enhancers_k562=read_data("enhancers_k562")

df_combined = pd.concat([
    df_promoters_hepg2.assign(dataset="promoters_hepg2"),
    df_promoters_k562.assign(dataset="promoters_k562"),
    df_enhancers_hepg2.assign(dataset="enhancers_hepg2"),
    df_enhancers_k562.assign(dataset="enhancers_k562")], ignore_index=True)

df_combined["cell_type"]=df_combined["dataset"].apply(lambda x: x.split("_")[1])
df_combined["re"]=df_combined["dataset"].apply(lambda x: x.split("_")[0])

# change track0/1 to cage
# change track2/3 to dhs
# change track4/5 to starr
# change track6/7 to sure
df_combined["track"]=df_combined["track"].apply(lambda x: "cage" if x=="track0" or x=="track1" else x)
df_combined["track"]=df_combined["track"].apply(lambda x: "dhs" if x=="track2" or x=="track3" else x)
df_combined["track"]=df_combined["track"].apply(lambda x: "starr" if x=="track4" or x=="track5" else x)
df_combined["track"]=df_combined["track"].apply(lambda x: "sure" if x=="track6" or x=="track7" else x)

# method 1:ks


track_list=[]
protein_list=[]
dstat_list=[]
for this_track in ["cage","dhs","starr","sure"]:
    df_this_track=df_combined[df_combined.track==this_track].reset_index(drop=True)
    for this_protein in df_this_track.protein.unique():
        logger.info(f"track:{this_track}, protein:{this_protein}")
        df_this_protein=df_this_track[df_this_track.protein==this_protein].reset_index(drop=True)
        if df_this_protein.shape[0]<10:
            continue
        dstat, _=ks_2samp(df_this_protein["ism_motif"],df_this_track["ism_motif"])
        # determine sign of dstat
        if df_this_protein["ism_motif"].median()<df_this_track["ism_motif"].median():
            dstat=-dstat
        track_list.append(this_track)
        protein_list.append(this_protein)
        dstat_list.append(dstat)
df_res=pd.DataFrame({"track":track_list,"protein":protein_list,"dstat":dstat_list})
df_res=df_res.pivot(index='protein', columns='track', values='dstat').reset_index()
df_res.to_csv("tf_ism_dstat_by_track_chip_true.csv")
            


# method 2: avg
# group by protein, track, avg the ism_motif
df_grouped=df_combined.groupby(["protein","track"]).agg({"ism_motif":"mean"}).reset_index()
df_pivot=df_grouped.pivot(index='protein', columns='track', values="ism_motif")
df_pivot.to_csv("tf_ism_avg_by_track_chip_true.csv")
#---------------------------------------------------------------------------
# Analysis2: analyze only the TFs presented in all datasets: compare between datasets
#---------------------------------------------------------------------------
# TODO: for each TF, for each track, find whether it is siginficantly emphasized in one dataset, one cell type, one RE

df=pd.read_csv("tf_individual_effect_by_file_chip_true.csv")
# pivot 
df_cage=df.pivot(index='protein', columns='dataset', values="avg_ism_cage")
# tfs unique to promoters_hepg2 are the row indices with NaN in promoters_k562, enhancers_hepg2, enhancers_k562
tfs_unique_promoters_hepg2=df_cage[df_cage["promoters_k562"].isnull() & df_cage["enhancers_hepg2"].isnull() & df_cage["enhancers_k562"].isnull()].index # IRF3
tfs_unique_promoters_k562=df_cage[df_cage["promoters_hepg2"].isnull() & df_cage["enhancers_hepg2"].isnull() & df_cage["enhancers_k562"].isnull()].index # SMAD2
tfs_unique_enhancers_hepg2=df_cage[df_cage["promoters_hepg2"].isnull() & df_cage["promoters_k562"].isnull() & df_cage["enhancers_k562"].isnull()].index # None
tfs_unique_enhancers_k562=df_cage[df_cage["promoters_hepg2"].isnull() & df_cage["promoters_k562"].isnull() & df_cage["enhancers_hepg2"].isnull()].index # HEY1
# tfs unique to enhancers are the row indices with NaN in promoters_k562, promoters_hepg2
tfs_unique_promoters=df_cage[df_cage["enhancers_hepg2"].isnull() & df_cage["enhancers_k562"].isnull()].index 
tfs_unique_enhancers=df_cage[df_cage["promoters_hepg2"].isnull() & df_cage["promoters_k562"].isnull()].index 
# remove nan rows
df_cage.dropna(inplace=True)

# nohup python3 get_tf_individual_effect.py > get_tf_individual_effect.out &








