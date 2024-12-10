import pandas as pd
from scipy.stats import ks_2samp
from loguru import logger


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from mutational_constraints import calc_constraint
from utils import get_track_num


# TODO: change ISM to ISA

prefix="/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_"


#------------------
# Helper functions
#------------------


def signed_ks_test(df,tf,isa_col):
    """
    calculate the ks test statistic for a given tf and a given importance
    """
    df_this_protein=df[df.protein==tf].reset_index(drop=True)
    if df_this_protein.shape[0]<10:
        return None
    dstat, _=ks_2samp(df_this_protein[isa_col],df[isa_col])
    # determine sign of dstat
    if df_this_protein[isa_col].median()<df[isa_col].median():
        dstat=-dstat
    return dstat


def calc_ks_stat_per_column(df,column):    
    protein_list=[]
    dstat_list=[]
    for this_protein in df.protein.unique():
        dstat=signed_ks_test(df,this_protein,column)
        protein_list.append(this_protein)
        dstat_list.append(dstat)
    df_res=pd.DataFrame({"protein":protein_list,f"dstat_{column}":dstat_list})
    return df_res





seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")



def process(df,cell_type):
    # calculate constraint
    df_constraints=calc_constraint(df,seq_extractor)
    # select and rename columns
    track_num_list=get_track_num(cell_type,classification=True)
    cols_retain=['protein']+[f'isa_track{i}' for i in track_num_list]
    df=df[cols_retain].copy()
    mapper={0:"cage_activity",
            1:"cage_activity",
            2:"dhs_activity",
            3:"dhs_activity",
            4:"starr_activity",
            5:"starr_activity",
            6:"sure_activity",
            7:"sure_activity",
            8:"cage_probability",
            9:"cage_probability",
            10:"dhs_probability",
            11:"dhs_probability",
            12:"starr_probability",
            13:"starr_probability",
            14:"sure_probability",
            15:"sure_probability"}
    df.rename(columns={f"isa_track{i}":f"isa_{mapper[i]}" for i in range(16)}, inplace=True)
    # calculate avg isa
    df_imp=df.groupby(["protein"]).mean().reset_index()
    df_imp.rename(columns={col:f"avg_{col}" for col in df_imp.columns[1:]}, inplace=True)
    # calculate ks
    for col in df.columns[1:]:
        df_res=calc_ks_stat_per_column(df,col)
        df_imp=pd.merge(df_imp,df_res,on="protein",how="outer")
    # merge with constraints
    df_imp=pd.merge(df_imp,df_constraints,on="protein",how="outer")
    return df_imp



def whole_analysis(df,col_name,dataset=None):
    if dataset is not None:
        df_sub=df[df[col_name]==dataset].reset_index(drop=True).copy()
    else:
        df_sub=df.copy()
    df_res=process(df_sub,dataset)
    df_res.to_csv(f"tf_effect_and_constraints_{dataset}.csv",index=False)


#---------------
# Analysis
#---------------
df_promoters_hepg2=pd.read_csv(f'{prefix}promoters_hepg2.csv')
df_promoters_k562=pd.read_csv(f'{prefix}promoters_k562.csv')
df_enhancers_hepg2=pd.read_csv(f'{prefix}enhancers_hepg2.csv')
df_enhancers_k562=pd.read_csv(f'{prefix}enhancers_k562.csv')

df = pd.concat([
    df_promoters_hepg2.assign(dataset="promoters_hepg2"),
    df_promoters_k562.assign(dataset="promoters_k562"),
    df_enhancers_hepg2.assign(dataset="enhancers_hepg2"),
    df_enhancers_k562.assign(dataset="enhancers_k562")
], ignore_index=True)

df["cell_type"]=df["dataset"].apply(lambda x: x.split("_")[1])

logger.info("Start analysis")

whole_analysis(df,"dataset")

whole_analysis(df,"dataset","promoters_hepg2")
logger.info("promoters_hepg2 done")
whole_analysis(df,"dataset","promoters_k562")
logger.info("promoters_k562 done")
whole_analysis(df,"dataset","enhancers_hepg2")
logger.info("enhancers_hepg2 done")
whole_analysis(df,"dataset","enhancers_k562")
logger.info("enhancers_k562 done")
whole_analysis(df,"cell_type","hepg2")
logger.info("hepg2 done")
whole_analysis(df,"cell_type","k562")
logger.info("k562 done")









# nohup python3 get_tf_individual_effect.py > get_tf_individual_effect.out &




