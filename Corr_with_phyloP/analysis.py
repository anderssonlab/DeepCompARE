import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import read_featimp,remove_nan_inf
from stat_tests import bin_and_label, calc_odds_ratio



# use ism by default

#----------------------
# Common functions
#----------------------

def get_track_num(file):
    if "k562" in file:
        return list(range(1,16,2))
    elif "hepg2" in file:
        return list(range(0,16,2))
    return None

#-------------------------------------------
# Count the number of conserved and accelerated
#-------------------------------------------


def read_ism_and_cons(file_prefix,phylop_prefix,track_num):
    # read in both files and merge into one dataframe
    df_ism=read_featimp(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/ism_{file_prefix}.csv",track_num)
    df_cons=pd.read_csv(f"{phylop_prefix}_{file_prefix}.csv",header=None)
    ism=df_ism.values.reshape(-1,1).squeeze()
    cons=df_cons.values.reshape(-1,1).squeeze()
    count_base=ism.shape[0]
    assert ism.shape==cons.shape
    ism,cons=remove_nan_inf(ism,cons)
    df = pd.DataFrame({'ism':ism,'cons':cons})
    # if cons<0, it is accelerated, else conserved
    df["evolution_type"]=df["cons"].apply(lambda x: "accelerated" if x<0 else "conserved")
    # if cons<-2.27 or > 2.27, it is significant, else non-significant
    df["sig"]=df["cons"].apply(lambda x: "Sig" if abs(x)>2.27 else "Non-sig")
    return df,count_base


def record_counts():
    dict_res={
        "file_prefix":[],
        "phylop":[],
        "count_base_total":[],
        "count_base_remain":[],
        "count_conserved":[],
        "count_accelerated":[],
        "count_sig_conserved":[],
        "count_sig_accelerated":[],
        "count_neg_ism":[],
        "count_pos_ism":[]
        }
    count_bases={}
    for file_prefix in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
        df_100way,count_bases["100way"]=read_ism_and_cons(file_prefix,"100way",0)
        df_241way,count_bases["241way"]=read_ism_and_cons(file_prefix,"241way",0)
        df_447way,count_bases["447way"]=read_ism_and_cons(file_prefix,"447way",0)
        df_447wayLRT,count_bases["447wayLRT"]=read_ism_and_cons(file_prefix,"447wayLRT",0)
        logger.info(f"Done with reading {file_prefix}")
        for phylop,df in zip(["100way","241way","447way","447wayLRT"],[df_100way,df_241way,df_447way,df_447wayLRT]):
            dict_res["file_prefix"].append(file_prefix)
            dict_res["phylop"].append(phylop)
            dict_res["count_base_total"].append(count_bases[phylop])
            dict_res["count_base_remain"].append(df.shape[0])
            dict_res["count_conserved"].append(df[df["evolution_type"]=="conserved"].shape[0])
            dict_res["count_accelerated"].append(df[df["evolution_type"]=="accelerated"].shape[0])
            dict_res["count_neg_ism"].append(df[df["ism"]<0].shape[0])
            dict_res["count_pos_ism"].append(df[df["ism"]>0].shape[0])
            df=df[df["sig"]=="Sig"].reset_index(drop=True)
            dict_res["count_sig_conserved"].append(df[df["evolution_type"]=="conserved"].shape[0])
            dict_res["count_sig_accelerated"].append(df[df["evolution_type"]=="accelerated"].shape[0])
    df_res=pd.DataFrame(dict_res)
    df_res.to_csv("counts2.csv",index=False)

    

if __name__=="__main__":
    record_counts()
    logger.info("Done with recording counts")


#-------------------------------------
# Plot odds ratio of significan bases
#-------------------------------------



# def read_ism_and_cons(file_prefix,phylop_prefix,track_num):
#     # read in both files and merge into one dataframe
#     df_ism=read_featimp(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/ism_{file_prefix}.csv",track_num)
#     df_cons=pd.read_csv(f"{phylop_prefix}_{file_prefix}.csv",header=None)
#     ism=df_ism.values.reshape(-1,1).squeeze()
#     cons=df_cons.values.reshape(-1,1).squeeze()
#     assert ism.shape==cons.shape
#     ism,cons=remove_nan_inf(ism,cons)
#     df = pd.DataFrame({'ism':ism,'cons':cons})
#     # split df to conserved and accelerated
#     df_conserved=df[df["cons"]>0].copy().reset_index(drop=True)
#     df_conserved=bin_and_label(df_conserved, 'cons', [0,2.27,np.inf])
#     df_conserved["Evolution_type"]=df_conserved.replace({"Bin":{"0-2.27":"Non-sig","2.27-inf":"Sig"}})["Bin"]
#     df_accelerated=df[df["cons"]<0].copy().reset_index(drop=True)
#     df_accelerated=bin_and_label(df_accelerated, 'cons', [-np.inf,-2.27,0])
#     df_accelerated["Evolution_type"]=df_accelerated.replace({"Bin":{"-inf--2.27":"Sig","-2.27-0":"Non-sig"}})["Bin"]
#     return df_conserved,df_accelerated


# def bin_ism_calc_OR(df,phylop_label,status_label):
#     df.drop(columns=["Bin"],inplace=True) 
#     df=bin_and_label(df, 'ism', [-np.inf,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,np.inf])
#     df_plot=df.groupby(["Bin","Evolution_type"]).size().unstack()
#     df_plot=df_plot[df_plot.sum(axis=1)>10]
#     df_plot[["or_Sig","pval"]]=df_plot.apply(lambda x: calc_odds_ratio(df_plot,x.name,"Sig"),axis=1).to_list()
#     df_plot["phyloP"]=phylop_label
#     df_plot["status"]=status_label
#     df_plot.drop(columns=["Non-sig","Sig"],inplace=True)
#     return df_plot



# def analyze_re(file_prefix,track_num):

#     df_100way_cons, df_100way_acc=read_ism_and_cons(file_prefix,"100way",track_num)
#     df_241way_cons, df_241way_acc=read_ism_and_cons(file_prefix,"241way",track_num)
#     df_447way_cons, df_447way_acc=read_ism_and_cons(file_prefix,"447way",track_num)
#     _, df_447wayLRT_acc=read_ism_and_cons(file_prefix,"447wayLRT",track_num)
    
#     or_100way_cons=bin_ism_calc_OR(df_100way_cons,"100way","conserved")
#     or_241way_cons=bin_ism_calc_OR(df_241way_cons,"241way","conserved")
#     or_447way_cons=bin_ism_calc_OR(df_447way_cons,"447way","conserved")

#     or_100way_acc=bin_ism_calc_OR(df_100way_acc,"100way","accelerated")
#     or_241way_acc=bin_ism_calc_OR(df_241way_acc,"241way","accelerated")
#     or_447way_acc=bin_ism_calc_OR(df_447way_acc,"447way","accelerated")
#     or_447wayLRT_acc=bin_ism_calc_OR(df_447wayLRT_acc,"447wayLRT","accelerated")
    
#     # concatenate the above by row
#     df_or=pd.concat([or_100way_cons,or_241way_cons,or_447way_cons,or_100way_acc,or_241way_acc,or_447way_acc,or_447wayLRT_acc],axis=0)
#     df_or["Predicted effect size"]=df_or.index
#     # negate or_Sig if status is accelerated
#     df_or["or_Sig"]=df_or.apply(lambda x: -x["or_Sig"] if x["status"]=="accelerated" else x["or_Sig"],axis=1)
#     df_or["alphas"]=df_or["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))
    
#     color_mapping = {'100way': "#1f77b4", '241way': '#ff7f0e', '447way': '#2ca02c', '447wayLRT': '#d62728'}

#     plt.figure(figsize=(6, 8))
#     for phylop, df_subset in df_or.groupby('phyloP'):
#         df_conserved=df_subset[df_subset["status"]=="conserved"]
#         df_accelerated=df_subset[df_subset["status"]=="accelerated"]
#         plt.plot(df_conserved['Predicted effect size'], df_conserved['or_Sig'], '--', color=color_mapping[phylop], label=phylop)
#         plt.scatter(df_conserved['Predicted effect size'], df_conserved['or_Sig'], color=color_mapping[phylop], alpha=df_subset['alphas'])
#         plt.plot(df_accelerated['Predicted effect size'], df_accelerated['or_Sig'], '--', color=color_mapping[phylop])
#         plt.scatter(df_accelerated['Predicted effect size'], df_accelerated['or_Sig'], color=color_mapping[phylop], alpha=df_subset['alphas'])
    
#     plt.xticks(rotation=45)
#     # put legend in the middle
#     plt.legend(loc='center right')  
#     # add horizontal line: y=0
#     plt.axhline(0, color='black', linestyle='--')
#     plt.ylabel("Odds ratio of significant base")
#     plt.title(f"Odds ratio in {file_prefix}, track {track_num}")
#     plt.savefig(f"OR_{file_prefix}_{track_num}.png")
#     plt.close()


# if __name__=="__main__":
#     logger.add("analysis.log")
#     for file_prefix in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
#         for track_num in get_track_num(file_prefix):
#             analyze_re(file_prefix,track_num)
#             logger.info(f"Done with {file_prefix}, track {track_num}")
    



# nohup python3 analysis.py &
