import pandas as pd
from scipy.stats import ks_2samp
from loguru import logger

prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info_and_motif_pair/motif_info_thresh_500_"


#------------------
# Heloer functions
#------------------
# TODO: remove everything about feat_imp

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
    # df_res=df_res[df_res["chip_evidence"]==True].reset_index(drop=True)
    return df_res


def calc_avg(df,column,new_col_prefix,dataset_annotation):
    """
    new_col_prefix can be either ism_motif or feat_imp_orig
    """
    df_res=df.groupby(["track","protein"])[column].mean().reset_index()
    df_res=df_res.pivot(index='protein', columns='track', values=column).reset_index()
    df_res.columns=["protein", f"avg_{new_col_prefix}_cage",f"avg_{new_col_prefix}_dhs",f"avg_{new_col_prefix}_starr",f"avg_{new_col_prefix}_sure"]
    df_res["dataset"]=dataset_annotation
    return df_res


def calc_ks_stat(df,column,new_col_prefix,dataset_annotation):
    track_list=[]
    protein_list=[]
    dstat_list=[]
    for this_track in df.track.unique():
        df_this_track=df[df.track==this_track].reset_index(drop=True)
        for this_protein in df_this_track.protein.unique():
            df_this_protein=df_this_track[df_this_track.protein==this_protein].reset_index(drop=True)
            if df_this_protein.shape[0]<10:
                continue
            dstat, _=ks_2samp(df_this_protein[column],df_this_track[column])
            # determine sign of dstat
            if df_this_protein[column].median()<df_this_track[column].median():
                dstat=-dstat
            track_list.append(this_track)
            protein_list.append(this_protein)
            dstat_list.append(dstat)
    df_res=pd.DataFrame({"track":track_list,"protein":protein_list,"dstat":dstat_list})
    df_res=df_res.pivot(index='protein', columns='track', values='dstat').reset_index()
    df_res.columns=["protein", f"dstat_{new_col_prefix}_cage",f"dstat_{new_col_prefix}_dhs",f"dstat_{new_col_prefix}_starr",f"dstat_{new_col_prefix}_sure"]
    df_res["dataset"]=dataset_annotation
    return df_res
              


def merge_multiple_files(df_list,on,how):
    df=df_list[0]
    for i in range(1,len(df_list)):
        df=pd.merge(df,df_list[i],on=on,how=how)
    return df


def calc_tf_importances(df,dataset):
    df_avg_ism = calc_avg(df, "ism_motif", "ism", dataset)
    df_dstat_ism = calc_ks_stat(df, "ism_motif", "ism", dataset)  
    return df_avg_ism,df_dstat_ism


def whole_analysis(df,subsetting_column,out_name):
    df_avg_ism_list = []
    df_dstat_ism_list = []
    for subsetting_value in df[subsetting_column].unique():
        logger.info(f"Processing {subsetting_value}")
        df_subset = df[df[subsetting_column]==subsetting_value].reset_index(drop=True)  
        df_avg_ism,df_dstat_ism = calc_tf_importances(df_subset,subsetting_value)
        df_avg_ism_list.append(df_avg_ism)
        df_dstat_ism_list.append(df_dstat_ism)
    df_avg_ism = pd.concat(df_avg_ism_list)
    df_dstat_ism = pd.concat(df_dstat_ism_list)
    df=merge_multiple_files([df_avg_ism,df_dstat_ism],
                on=["protein","dataset"],
                how="outer")
    df.to_csv(out_name,index=False)



#---------------
# Analysis
#---------------

df_promoters_hepg2=read_data("promoters_hepg2")
logger.info("promoters_hepg2 read")
df_promoters_k562=read_data("promoters_k562")
logger.info("promoters_k562 read")
df_enhancers_hepg2=read_data("enhancers_hepg2")
logger.info("enhancers_hepg2 read")
df_enhancers_k562=read_data("enhancers_k562")
logger.info("enhancers_k562 read")

df_combined = pd.concat([
    df_promoters_hepg2.assign(dataset="promoters_hepg2"),
    df_promoters_k562.assign(dataset="promoters_k562"),
    df_enhancers_hepg2.assign(dataset="enhancers_hepg2"),
    df_enhancers_k562.assign(dataset="enhancers_k562")
], ignore_index=True)

df_combined["cell_type"]=df_combined["dataset"].apply(lambda x: x.split("_")[1])
df_combined["re"]=df_combined["dataset"].apply(lambda x: x.split("_")[0])

whole_analysis(df_combined,"dataset","tf_individual_effect_by_file.csv")
logger.info("aggregation by file done")
whole_analysis(df_combined,"cell_type","tf_individual_effect_by_cell_type.csv")
logger.info("aggregation by cell type done")
# whole_analysis(df_combined,"re","tf_individual_effect_by_re.csv")
# logger.info("aggregation by re done")







# nohup python3 get_tf_individual_effect.py > get_tf_individual_effect.out &








#----------------------------------------------------------------
# Optional: test correlation between ism_motif and feat_imp_orig for count_TF=1 and >1
# Conclusion: No difference
# TODO: does count_tf>1 cause a systematic drop in ism_motif and feat_imp_orig?
#----------------------------------------------------------------
# def test_corr(df):
#     ism_nonredundant=df[df.count_TF_no_thresh==1]["ism_motif"]
#     ism_redundant=df[df.count_TF_no_thresh>1]["ism_motif"]
#     gradxinp_nonredundant=df[df.count_TF_no_thresh==1]["feat_imp_orig"]
#     gradxinp_redundant=df[df.count_TF_no_thresh>1]["feat_imp_orig"]
#     print("nonredundant",ism_nonredundant.corr(gradxinp_nonredundant))
#     print("redundant",ism_redundant.corr(gradxinp_redundant))

# test_corr(df_promoters_hepg2)
# test_corr(df_promoters_k562)
# test_corr(df_enhancers_hepg2)
# test_corr(df_enhancers_k562)
#----------------------------------------------------------------


