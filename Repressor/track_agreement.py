import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import get_track_num
from mutational_constraints import calc_constraint
from seq_ops import SeqExtractor
# ---------------
# Helper functions
# ---------------

def read_file(file_suffix):
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_suffix}.csv")
    df=df[df['chip_evidence']==True].reset_index(drop=True)
    track_num_list=get_track_num(file_suffix,classification=False)
    cols_retain=['chromosome', 'start', 'end', 'protein', '241way','af']+\
        [f"ism_track{track_num}" for track_num in track_num_list]
    df=df[cols_retain].copy()
    df.drop_duplicates(subset=["chromosome","start","end","protein"], inplace=True)
    mapper={0:"cage",
            1:"cage",
            2:"dhs",
            3:"dhs",
            4:"starr",
            5:"starr",
            6:"sure",
            7:"sure"}
    df.rename(columns={f"ism_track{i}":f"ism_{mapper[i]}" for i in range(8)}, inplace=True)
    df["max_af"] = df["af"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["min_af"] = df["af"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["num_variants"] = df["af"].apply(lambda x: len(str(x).split(":")))
    df["num_common_variants"] = df["af"].apply(lambda x: sum([float(i)>0.001 for i in str(x).split(":")]))
    df["num_rare_variants"] = df["af"].apply(lambda x: sum([float(i)<=0.001 for i in str(x).split(":")]))
    df["241way_max"] = df["241way"].apply(lambda x: max([float(i) for i in str(x).split(":")]))
    df["241way_min"] = df["241way"].apply(lambda x: min([float(i) for i in str(x).split(":")]))
    df["dataset"]=file_suffix
    return df


def detect_conflict(df,threshold):
    """
    Add 1 column called "conflicting" to the dataframe
    remove rows with no conclusions
    """
    df["repressor_cage"]=df["ism_cage"]< -threshold
    df["activator_cage"]=df["ism_cage"]> threshold
    df["repressor_dhs"]=df["ism_dhs"]< -threshold
    df["activator_dhs"]=df["ism_dhs"]> threshold
    df["repressor_starr"]=df["ism_starr"]< -threshold
    df["activator_starr"]=df["ism_starr"]> threshold
    df["repressor_sure"]=df["ism_sure"]< -threshold
    df["activator_sure"]=df["ism_sure"]> threshold
    df["num_repressor"]=df[["repressor_cage","repressor_dhs","repressor_starr","repressor_sure"]].sum(axis=1)
    df["num_activator"]=df[["activator_cage","activator_dhs","activator_starr","activator_sure"]].sum(axis=1)
    # remove rows without conclusions
    df=df[(df["num_repressor"]>0) | (df["num_activator"]>0)].reset_index(drop=True)
    # remove rows with conflictions
    df["conflicting"] = (df["num_repressor"] > 0) & (df["num_activator"] > 0)
    df["activator"]=df["num_activator"]>0
    df["repressor"]=df["num_repressor"]>0
    # drop the columns
    df.drop(columns=["repressor_cage","activator_cage",
                     "repressor_dhs","activator_dhs",
                     "repressor_starr","activator_starr",
                     "repressor_sure","activator_sure",
                     "num_repressor","num_activator"], inplace=True)
    return df
    

seq_extractor = SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")


# ---------------
# read data
# ---------------
df_promoters_hepg2=read_file("promoters_hepg2")
df_promoters_k562=read_file("promoters_k562")
df_enhancers_hepg2=read_file("enhancers_hepg2")
df_enhancers_k562=read_file("enhancers_k562")
# concat, and assign new labels in the column "dataset"
df=pd.concat([df_promoters_hepg2,df_promoters_k562,df_enhancers_hepg2,df_enhancers_k562],
             axis=0)
df=detect_conflict(df,0.05)
df=df[df["conflicting"]==False].reset_index(drop=True)




#--------------------
# task1: group by protein and dataset, count percentage of repressors, and total occurences
#--------------------
df["repressor_count"]=df["repressor"].astype(int)
df["activator_count"]=df["activator"].astype(int)
df_fraction_repressor=df.groupby(["protein","dataset"]).agg({"repressor":"mean",
                                                 "conflicting":"count",
                                                 "repressor_count":"sum",
                                                 "activator_count":"sum"
                                                 }).reset_index()
df_fraction_repressor.rename(columns={"conflicting":"num_samples",
                     "repressor":"fraction_repressor"
                   }, inplace=True)
df_fraction_repressor=df_fraction_repressor[df_fraction_repressor["num_samples"]>10].reset_index(drop=True)
df_fraction_repressor["min_count"]=df_fraction_repressor[["activator_count","repressor_count"]].min(axis=1)

sns.histplot(data=df_fraction_repressor, x="fraction_repressor",bins=50)
plt.savefig("hist_fraction_repressor.pdf")
plt.close()

df_fraction_repressor[df_fraction_repressor["fraction_repressor"]>0.7].reset_index(drop=True)
df_fraction_repressor[df_fraction_repressor["fraction_repressor"]==0].reset_index(drop=True)

df_fraction_repressor[df_fraction_repressor["dataset"]=="enhancers_hepg2"].sort_values(by="min_count",ascending=False).reset_index(drop=True)

#--------------------------------------------------
# task2: for each tf, compare their constraints  between repressor and activator state
#--------------------------------------------------
def calculate_constraint_by_state(df,seq_extractor):
    df_constraint=calc_constraint(df,seq_extractor)
    df_constraint["state"]=df_constraint.protein.str.split("_").str[-1]
    df_constraint["protein"]=df_constraint.protein.str.split("_").str[0]
    df_constraint=df_constraint.pivot(index="protein", columns="state", values="z").reset_index()
    return df_constraint


# if activator, then concat "_activator" to the protein name
df["protein"]=df.apply(lambda x: x["protein"]+"_activator" if x["activator"] else x["protein"], axis=1)
df["protein"]=df.apply(lambda x: x["protein"]+"_repressor" if x["repressor"] else x["protein"], axis=1)




df_constraint_promoters_hepg2=calculate_constraint_by_state(df[df["dataset"]=="promoters_hepg2"],seq_extractor)
df_constraint_promoters_hepg2["dataset"]="promoters_hepg2"
df_constraint_promoters_k562=calculate_constraint_by_state(df[df["dataset"]=="promoters_k562"],seq_extractor)
df_constraint_promoters_k562["dataset"]="promoters_k562"
df_constraint_enhancers_hepg2=calculate_constraint_by_state(df[df["dataset"]=="enhancers_hepg2"],seq_extractor)
df_constraint_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_constraint_enhancers_k562=calculate_constraint_by_state(df[df["dataset"]=="enhancers_k562"],seq_extractor)
df_constraint_enhancers_k562["dataset"]="enhancers_k562"


df_constraint=pd.concat([df_constraint_promoters_hepg2,df_constraint_promoters_k562,df_constraint_enhancers_hepg2,df_constraint_enhancers_k562],axis=0)


# merge with df_fraction_repressor, by protein and dataset
df_constraint=pd.merge(df_constraint, df_fraction_repressor, on=["protein","dataset"], how="inner")
# plot activator vs repressor, hue=fraction_repressor, size=min_count
for file in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
    df_subset=df_constraint[df_constraint["dataset"]==file].reset_index(drop=True)
    plt.figure(figsize=(7,7))
    sns.scatterplot(data=df_subset,x="activator", y="repressor",hue="fraction_repressor", size="min_count")
    # annotate the points
    for i in range(df_subset.shape[0]):
        if df_subset["min_count"][i]>44 or df_subset["activator"][i]>5 or df_subset["repressor"][i]<-3:
            plt.text(df_subset["activator"][i], df_subset["repressor"][i], df_subset["protein"][i], fontsize=6)
    # plot diagonal
    min_val=min(df_subset["activator"].min(),df_subset["repressor"].min())
    max_val=max(df_subset["activator"].max(),df_subset["repressor"].max())
    plt.plot([min_val,max_val],[min_val,max_val], color="red")
    plt.title(f"{file}")
    plt.savefig(f"scatter_constraint_activator_vs_repressor_{file}.pdf")
    plt.close()



#--------------------------------------------------
# task3: for each tf, compare their number of variants/evolutionary conservation between repressor and activator state
#--------------------------------------------------
df_grouped=df.groupby(["protein","dataset"]).agg({"num_rare_variants":"mean",
                                                  "241way_max":"mean",
                                                  "num_common_variants":"mean"
                                                  }).reset_index()
df_grouped["state"]=df_grouped.protein.str.split("_").str[-1]
df_grouped["protein"]=df_grouped.protein.str.split("_").str[0]


def pivot_and_plot(df_orig,file_suffix,var):
    df=df_orig.copy()
    df=df[df["dataset"]==file_suffix].reset_index(drop=True)
    df=df.pivot(index="protein", columns="state", values=var).reset_index()
    df_fraction=df_fraction_repressor[df_fraction_repressor['dataset']==file_suffix].copy()
    # merge with df_fraction_repressor
    df=pd.merge(df, df_fraction, on="protein", how="inner")
    plt.figure(figsize=(7,7))
    sns.scatterplot(data=df, x="activator", y="repressor", hue="fraction_repressor", size="min_count")
    # text annotation
    texts=[]
    for i in range(df.shape[0]):
        if df["min_count"][i]>44 or abs(df["activator"][i]-df["repressor"][i])>0.2:
            texts.append(plt.text(df["activator"][i], df["repressor"][i], df["protein"][i], fontsize=6))
    # plot diagonal
    adjust_text(texts)
    min_val=min(df["activator"].min(),df["repressor"].min())
    max_val=max(df["activator"].max(),df["repressor"].max())
    plt.plot([min_val,max_val],[min_val,max_val], color="red")
    plt.title(f"{var} of {file_suffix}")
    plt.savefig(f"scatter_{var}_activator_vs_repressor_{file_suffix}.pdf")
    plt.close()



pivot_and_plot(df_grouped,"promoters_hepg2","num_common_variants")
pivot_and_plot(df_grouped,"promoters_k562","num_common_variants")
pivot_and_plot(df_grouped,"enhancers_hepg2","num_common_variants")
pivot_and_plot(df_grouped,"enhancers_k562","num_common_variants")



pivot_and_plot(df_grouped,"promoters_hepg2","241way_max")
pivot_and_plot(df_grouped,"promoters_k562","241way_max")
pivot_and_plot(df_grouped,"enhancers_hepg2","241way_max")
pivot_and_plot(df_grouped,"enhancers_k562","241way_max")



# ---------------
# Archived
# ---------------


# # distribution of isa scores
# for file_suffix in ["promoters_hepg2","promoters_k562","enhancers_hepg2","enhancers_k562"]:
#     df_sub=df[df["dataset"]==file_suffix].reset_index(drop=True).copy()
#     df_long=pd.melt(df_sub, id_vars=["protein"], value_vars=[f"ism_{track_name}" for track_name in ["cage","dhs","starr","sure"]], var_name="track", value_name="ism_score")
#     sns.histplot(data=df_long, x="ism_score", hue="track", bins=100)
#     # y log scale
#     plt.yscale("log")
#     plt.title(f"ISA score distribution for {file_suffix}")
#     plt.savefig(f"hist_isa_distribution_{file_suffix}.pdf")
#     plt.close()





# # do tracks agree on repressor/activator status?
# def threshold_vs_conflict(df,threshold):
#     df=detect_conflict(df,threshold)
#     return df["conflicting"].sum()/df.shape[0],df.shape[0]


# thresh_list=np.linspace(0,0.2,200)
# res=[threshold_vs_conflict(df, thresh) for thresh in thresh_list]
# # split the tuple
# conflict_list=[x[0] for x in res]
# num_samples=[x[1] for x in res]
# plt.plot(thresh_list, conflict_list)
# plt.xlabel("Threshold")
# plt.ylabel("Fraction of conflicting predictions")
# plt.savefig("threshold_vs_conflict.pdf")
# plt.close()


# plt.plot(thresh_list, num_samples,label="Number of samples")
# plt.xlabel("Threshold")
# plt.ylabel("Number of samples")
# plt.savefig("threshold_vs_num_samples.pdf")
# plt.close()

