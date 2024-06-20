import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from scipy.stats import ks_2samp
from matplotlib.colors import LinearSegmentedColormap

# Define a custom color map
colors = ["blue", "lightgrey", "yellow"]  # blue for negative, lightgray for zero, red for positive
n_bins = 100  # Increase this number to make the transition smoother
cmap_name = "custom"
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)



#---------------------------------------
# Functions for downstream analysis
#---------------------------------------
def get_condition_percentage(group, condition, colname):
    condition_rows = condition(group)
    return pd.Series({colname: condition_rows})


def read_file(file_name):
    df=pd.read_csv(file_name)
    df=df[(df["ism_score_mut1"]>0) & (df["ism_score_mut2"]>0)].reset_index(drop=True)
    df["ism2_wo_protein1"]=df['pred_mut1']-df['pred_mut_both']
    df["ism2_wo_protein1_null"]=df['pred_mut1_null']-df['pred_mut_both_null']
    df["diff"]=df["ism_score_mut2"]-df["ism2_wo_protein1"]
    df["diff_null"]=df["ism_score_mut2"]-df["ism2_wo_protein1_null"]
    df["distance"]=np.abs(df["start_rel2"]-df["start_rel1"])
    df["distance_null"]=np.abs(df["start_rel2"]-df["start_rel_null"])
    return df


def get_track_num(dataset):
    if "hepg2" in dataset:
        return [0,2,4,6]
    elif "k562" in dataset:
        return [1,3,5,7]
    else:
        return None
#---------------------------------------
# analysis 1: get df_null_of_balance
#---------------------------------------


def calculate_threshold(df,colname,thresh):
    # assume index is distance
    # smooth curve
    df=df.copy()
    df["value_smoothed"]=df[colname].rolling(window=5).mean()
    df_sub=df[df[colname]<thresh]
    if df_sub.shape[0]>0:
        return df_sub.index[0]
    else:
        return df.index[-1]

def get_threshold_distance(df_input,diff_thresh=1e-4):
    df=df_input.copy()
    df["diff_abs"]=np.abs(df["diff"])
    df["diff_null_abs"]=np.abs(df["diff_null"])
    df_distance=df.groupby(["distance"]).agg({"diff_abs":"mean"})
    df_distance_null=df.groupby(["distance_null"]).agg({"diff_null_abs":"mean"})
    # find the distance where smooth diff drop below 1e-4
    dist_thresh1=calculate_threshold(df_distance,"diff_abs",diff_thresh)
    dist_thresh2=calculate_threshold(df_distance_null,"diff_null_abs",diff_thresh)
    return max(dist_thresh1,dist_thresh2)

def get_additivity_threshold(df):
    # get 99 percentile of posivite values in diff_null as sup threshold
    # get 1 percentile of negative values in diff_null as sub threshold
    if df["diff_null"].max()<=0:
        super_thersh=0
    else:    
        diff_null_pos=df[df["diff_null"]>0]["diff_null"]
        super_thersh=np.percentile(diff_null_pos,95)
    if df["diff_null"].min()>=0:
        sub_thersh=0
    else:
        diff_null_neg=df[df["diff_null"]<0]["diff_null"]
        sub_thersh=np.percentile(diff_null_neg,5)
    return sub_thersh,super_thersh 

def ks_test(df):
    df_super=df[df["super-additivity"]]
    df_sub=df[~df["super-additivity"]]
    d,p=ks_2samp(df_super["distance"],df_sub["distance"])
    sign=df_sub["distance"].median()-df_super["distance"].median()
    d=d*sign
    return d,p,df_sub["distance"].median(),df_super["distance"].median()

    
def analyze_one_file(file_name):
    df=read_file(file_name)
    
    # lists to save info
    tf_list=[]
    wilcoxon_null_list=[]
    wilcoxon_list=[]
    median_null_list=[]
    median_list=[]
    threshold_distance_list=[]
    sub_thresh_list=[]
    super_thresh_list=[]
    total_count_list=[] # record the number of counts within distance threshold
    total_pos_count_list=[] # record the number of counts within distance threshold and diff>0
    total_neg_count_list=[] # record the number of counts within distance threshold and diff<0
    sig_super_count_list=[] # record the number of counts with significant superadditivity
    sig_sub_count_list=[] # record the number of counts with significant subadditivity
    super_sub_ratio_list=[] # record the ratio of significant superadditivity to significant subadditivity
    
    dstat_list=[]
    pval_list=[]
    distance_median_sub_list=[]
    distance_median_super_list=[]
    
    for protein in df["protein2"].unique():
        
        df_sub=df[df["protein2"]==protein].reset_index(drop=True)
        if df_sub.shape[0]<10:
            continue
        threshold_distance=get_threshold_distance(df_sub)
        df_sub=df_sub[df_sub["distance"]<threshold_distance].reset_index(drop=True)
        if df_sub.shape[0]<10:
            continue
        
        df_sub["super-additivity"]=(df_sub["diff"]>0)
        if df_sub["super-additivity"].nunique()==1:
            continue
        
        # wilcoxon to test whether median is significantly nonzero
        _,p=wilcoxon(df_sub["diff"])
        _,p_null=wilcoxon(df_sub["diff_null"])
        
        # get sub/super additivity threshold 
        sub_thresh,super_thresh=get_additivity_threshold(df_sub)
        sig_super_occurances=df_sub[(df_sub["diff"]>super_thresh)].shape[0]
        sig_sub_occurances=df_sub[(df_sub["diff"]<sub_thresh)].shape[0]
        super_sub_ratio=np.log2((sig_super_occurances+1)/(sig_sub_occurances+1)) # pseudocount
        
        # KS test for distance
        dstat,pval,distance_median_sub,distance_median_super=ks_test(df_sub)
        dstat_list.append(dstat)
        pval_list.append(pval)
        distance_median_sub_list.append(distance_median_sub)
        distance_median_super_list.append(distance_median_super)
        
        # add to list
        tf_list.append(protein)
        wilcoxon_null_list.append(p_null)
        wilcoxon_list.append(p)
        median_null_list.append(df_sub["diff_null"].median())
        median_list.append(df_sub["diff"].median())
        threshold_distance_list.append(threshold_distance)
        sub_thresh_list.append(sub_thresh)
        super_thresh_list.append(super_thresh)
        total_count_list.append(df_sub.shape[0])
        total_pos_count_list.append(df_sub[df_sub["diff"]>0].shape[0])
        total_neg_count_list.append(df_sub[df_sub["diff"]<0].shape[0])
        sig_super_count_list.append(sig_super_occurances)
        sig_sub_count_list.append(sig_sub_occurances)
        super_sub_ratio_list.append(super_sub_ratio)

    df_res=pd.DataFrame({"protein":tf_list,
                              "median":median_list,
                              "p":wilcoxon_list,
                              "median_null":median_null_list,
                              "p_null":wilcoxon_null_list,
                              "threshold_distance":threshold_distance_list,
                              "sub_thresh":sub_thresh_list,
                              "super_thresh":super_thresh_list,
                              "total_count":total_count_list,
                              "total_pos_count":total_pos_count_list,
                              "total_neg_count":total_neg_count_list,
                              "sig_super_count":sig_super_count_list,
                              "sig_sub_count":sig_sub_count_list,
                              "super_sub_ratio":super_sub_ratio_list,
                              "dstat":dstat_list,
                              "pval":pval_list,
                              "distance_median_sub":distance_median_sub_list,
                              "distance_median_super":distance_median_super_list})
    return df_res

    
# analysis for all datasets
df_list=[]
for dataset in ["enhancers_hepg2","promoters_hepg2"]:
    for track_num in get_track_num(dataset):
        df=analyze_one_file(f"Pd1_df_mutate_pair/df_mutate_pair_{dataset}_remap_Hep-G2,track{track_num}.csv")
        df["dataset"]=dataset
        df["track_num"]=track_num
        df_list.append(df)
        
for dataset in ["enhancers_k562","promoters_k562"]:
    for track_num in get_track_num(dataset):
        df=analyze_one_file(f"Pd1_df_mutate_pair/df_mutate_pair_{dataset}_remap_K-562,track{track_num}.csv")
        df["dataset"]=dataset
        df["track_num"]=track_num
        df_list.append(df)
        
df_all=pd.concat(df_list)
df_all.to_csv("df_sup_sub.csv" ,index=False)    






#---------------------------------------
# analysis 2: plot wilcoxon (median) test results for null and non-null
# median < 0: ism_score_mut2 < ism2_wo_protein1: subadditivity
#---------------------------------------
df=pd.read_csv("df_sup_sub.csv")

(df["median"]<0).sum() # 900/2637
(df["median_null"]<0).sum() # 1441/2637
(df.p>0.05).sum() # 583
(df.p_null>0.05).sum() # 1112

df_sig_sup_null = df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: ((group['median_null'] > 0) & (group['p_null'] < 0.05)).mean(),colname="sig_super_null").reset_index()
df_sig_sub_null = df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: ((group['median_null'] < 0) & (group['p_null'] < 0.05)).mean(),colname="sig_sub_null").reset_index()
df_sig_sup=df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: ((group['median'] > 0) & (group['p'] < 0.05)).mean(),colname="sig_super").reset_index()
df_sig_sub=df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: ((group['median'] < 0) & (group['p'] < 0.05)).mean(),colname="sig_sub").reset_index()

# merge 4 dfs by dataset and track number
merged_df = pd.merge(df_sig_sub, df_sig_sup, on=['dataset', 'track_num'], how='inner')
merged_df = pd.merge(merged_df, df_sig_sub_null, on=['dataset', 'track_num'], how='inner')
merged_df = pd.merge(merged_df, df_sig_sup_null, on=['dataset', 'track_num'], how='inner')

merged_df['track_num'] = merged_df['track_num'].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})

# is any track biased?
plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_sub', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig('Plots/wilcoxon_sig_sub.pdf')
plt.close()

plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_sub_null', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig('Plots/wilcoxon_sig_sub_null.pdf')
plt.close()

plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_super', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig('Plots/wilcoxon_sig_super.pdf')
plt.close()

plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_super_null', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig('Plots/wilcoxon_sig_super_null.pdf')
plt.close()


#-----------------------------------------------------------------------------
# analysis 3: plot significant sub/super additivity
# reproduce file (x)-sub/superpercentage(y)-track(hue) plot using super_sub_ratio
#-----------------------------------------------------------------------------

df=pd.read_csv("df_sup_sub.csv")
df=df[(df["sig_super_count"]+df["sig_sub_count"])>10].reset_index(drop=True)

(df.super_sub_ratio<0).mean() # 0.37

threshold=0

df_sig_sup=df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: (group['super_sub_ratio'] > threshold).mean(),colname="sig_super").reset_index()
df_sig_sub=df.groupby(['dataset', 'track_num']).apply(get_condition_percentage, condition=lambda group: (group['super_sub_ratio'] < -threshold).mean(),colname="sig_sub").reset_index()

merged_df = pd.merge(df_sig_sub, df_sig_sup, on=['dataset', 'track_num'], how='inner')
merged_df['track_num'] = merged_df['track_num'].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})

plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_sub', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig(f'Plots/thresholded_sig_sub_{threshold}.pdf')
plt.close()

plt.figure(figsize=(8,6))
sns.catplot(x='dataset', y='sig_super', hue='track_num', data=merged_df)
# rotate x labels by 45 degrees
plt.xticks(rotation=45)
# make lower margin larger
plt.subplots_adjust(bottom=0.28)
plt.savefig(f'Plots/thresholded_sig_super_{threshold}.pdf')
plt.close()


#-----------------------------------------------------
# Analysis 4: do tracks show difference in threshold distance
#-----------------------------------------------------

df=pd.read_csv("df_sup_sub.csv")
df['track_num'] = df['track_num'].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})

# seaborn violin plot to show distribuion of threshold_distance
sns.boxplot(x='dataset', y='threshold_distance', hue='track_num', data=df)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.28)
plt.savefig('Plots/distribution_of_threshold_distance.pdf')
plt.close()




#----------------------------------------------------------------------------------------------------
# Analysis 5: do superadditivity require closer distance than subadditivity?
#----------------------------------------------------------------------------------------------------

    
df=pd.read_csv("df_sup_sub.csv")
df=df[(df["sig_super_count"]+df["sig_sub_count"])>10].reset_index(drop=True)
df["track_num"]=df["track_num"].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})


(df.dstat>0).mean() # 82%, ignoring dataset and track_num
# how many rows have dstat>0 and pval<0.05
((df.dstat>0) & (df.pval<0.05)).mean()# 75%, ignoring dataset and track_num

df_distance_sub=df[["protein","dataset","track_num","distance_median_sub"]].copy().rename(columns={"distance_median_sub":"distance_median"})
df_distance_sub["additivity"]="sub-additive"
df_distance_super=df[["protein","dataset","track_num","distance_median_super"]].copy().rename(columns={"distance_median_super":"distance_median"})
df_distance_super["additivity"]="super-additive"
df_distance=pd.concat([df_distance_sub,df_distance_super])

# violin plot to show distribution of distance_median, but lose information about pairing 
for assay in ["cage","dhs","starr","sure"]:
    # violin plot to show distribution of distance_median, but lose information about pairing 
    sns.violinplot(x="dataset",y="distance_median",hue="additivity",data=df_distance[df_distance.track_num==assay],split=True,inner="quart")
    plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.28)
    plt.title(f"Cooperitivity measured by {assay}")
    plt.savefig(f'Plots_sup_sub_distance/distance_median_{assay}.pdf')
    plt.close()
    sns.scatterplot(x="distance_median_sub",y="distance_median_super",hue="dataset",data=df[df.track_num==assay])
    plt.title(f"Cooperitivity measured by {assay}")
    min_val=min(df["distance_median_sub"].min(),df["distance_median_super"].min())
    max_val=max(df["distance_median_sub"].max(),df["distance_median_super"].max())
    plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
    plt.savefig(f'Plots_sup_sub_distance/distance_pair_{assay}.pdf')
    plt.close()


#----------------------------------------------------------------------------------------------------
# Analysis 6: plot heatmap of log2(#sup/#sub) 
#----------------------------------------------------------------------------------------------------
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform

df=pd.read_csv("df_sup_sub.csv")
df=df[(df["sig_super_count"]+df["sig_sub_count"])>10].reset_index(drop=True)
df=df[df.sig_sub_count>10].reset_index(drop=True)
df=df[df.sig_super_count>10].reset_index(drop=True)
df["track_num"]=df["track_num"].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})


for assay in ["cage","dhs","starr","sure"]:
    df_ratio=df[df.track_num==assay].pivot(index="dataset", columns="protein", values="super_sub_ratio")
    nan_idx=np.isnan(df_ratio)
    df_ratio = df_ratio.fillna(0)
    # Perform hierarchical clustering on the filled data
    row_linkage = linkage(squareform(pdist(df_ratio)), method='average')
    col_linkage = linkage(squareform(pdist(df_ratio.T)), method='average')
    # Get the leaves order
    df_ratio = df_ratio.iloc[leaves_list(row_linkage), :]
    df_ratio = df_ratio.iloc[:, leaves_list(col_linkage)]
    plt.figure(figsize=(50, 3))
    sns.heatmap(df_ratio, cmap=custom_cmap, cbar_kws={'label': 'log2(#sup/#sub)'},xticklabels=True)
    plt.subplots_adjust(bottom=0.5)
    plt.title(f"Heatmap of log2(#sup/#sub) measured by {assay}")
    plt.savefig(f"Plots/heatmap_{assay}.pdf",dpi=300)
    plt.close()


# nohup python3 get_df_sub_super.py > get_df_sub_super.out &