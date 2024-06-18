import pandas as pd


prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_"

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
    df_res=pd.concat(df_list)
    # df_res=df_res[df_res["chip_evidence"]==True].reset_index(drop=True)
    return df_res

def calc_avg_ism(df):
    df_res=df.groupby(["track","protein"])["ism_motif"].mean().reset_index()
    df_res=df_res.pivot(index='protein', columns='track', values='ism_motif').reset_index()
    df_res.columns=["protein", "ism_cage", "ism_dhs", "ism_starr", "ism_sure"]
    return df_res

def calc_avg_gradxinp(df):
    df_res=df.groupby(["track","protein"])["feat_imp_orig"].mean().reset_index()
    df_res=df_res.pivot(index='protein', columns='track', values='feat_imp_orig').reset_index()
    df_res.columns=["protein","gradxinp_cage","gradxinp_dhs","gradxinp_starr","gradxinp_sure"]
    return df_res



df_promoters_hepg2=read_data("promoters_hepg2")
df_promoters_k562=read_data("promoters_k562")
df_enhancers_hepg2=read_data("enhancers_hepg2")
df_enhancers_k562=read_data("enhancers_k562")



#----------------------------------------------------------------
# Optional: test correlation between ism_motif and feat_imp_orig for count_TF=1 and >1
# Conclusion: No difference
#----------------------------------------------------------------
def test_corr(df):
    ism_nonredundant=df[df.count_TF_no_thresh==1]["ism_motif"]
    ism_redundant=df[df.count_TF_no_thresh>1]["ism_motif"]
    gradxinp_nonredundant=df[df.count_TF_no_thresh==1]["feat_imp_orig"]
    gradxinp_redundant=df[df.count_TF_no_thresh>1]["feat_imp_orig"]
    print("nonredundant",ism_nonredundant.corr(gradxinp_nonredundant))
    print("redundant",ism_redundant.corr(gradxinp_redundant))

test_corr(df_promoters_hepg2)
test_corr(df_promoters_k562)
test_corr(df_enhancers_hepg2)
test_corr(df_enhancers_k562)
#----------------------------------------------------------------





df_ism_promoters_hepg2=calc_avg_ism(df_promoters_hepg2)
df_ism_promoters_k562=calc_avg_ism(df_promoters_k562)
df_ism_enhancers_hepg2=calc_avg_ism(df_enhancers_hepg2)
df_ism_enhancers_k562=calc_avg_ism(df_enhancers_k562)

df_ism_promoters_hepg2["dataset"]="promoters_hepg2"
df_ism_promoters_k562["dataset"]="promoters_k562"
df_ism_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_ism_enhancers_k562["dataset"]="enhancers_k562"


df_gradxinp_promoters_hepg2=calc_avg_gradxinp(df_promoters_hepg2)
df_gradxinp_promoters_k562=calc_avg_gradxinp(df_promoters_k562)
df_gradxinp_enhancers_hepg2=calc_avg_gradxinp(df_enhancers_hepg2)
df_gradxinp_enhancers_k562=calc_avg_gradxinp(df_enhancers_k562)

df_gradxinp_promoters_hepg2["dataset"]="promoters_hepg2"
df_gradxinp_promoters_k562["dataset"]="promoters_k562"
df_gradxinp_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_gradxinp_enhancers_k562["dataset"]="enhancers_k562"


# merge all datasets by protein and dataset
df_ism=pd.concat([df_ism_promoters_hepg2,df_ism_promoters_k562,df_ism_enhancers_hepg2,df_ism_enhancers_k562])
df_gradxinp=pd.concat([df_gradxinp_promoters_hepg2,df_gradxinp_promoters_k562,df_gradxinp_enhancers_hepg2,df_gradxinp_enhancers_k562])
df=pd.merge(df_ism,df_gradxinp,on=["protein","dataset"],how="inner")
df.to_csv("tf_individual_effect.csv",index=False)