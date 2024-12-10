import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,spearmanr
from adjustText import adjust_text
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import bin_and_label
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_index, calculate_tf_cooperativity_index


# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------

def read_all_files(threshold):
    # read cooperativity
    df_promoter_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_promoters_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_enhancer_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_hepg2.csv",threshold=threshold,track_nums=[0,2,4,6])
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_enhancers_k562.csv",threshold=threshold,track_nums=[1,3,5,7])
    df_promoter_hepg2["dataset"]="promoter_hepg2"
    df_promoter_k562["dataset"]="promoter_k562"
    df_enhancer_hepg2["dataset"]="enhancer_hepg2"
    df_enhancer_k562["dataset"]="enhancer_k562"
    df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562],axis=0)
    return df


def write_ci_by_cell_type(df_orig,cell_type):
    df=df_orig.copy()
    if cell_type!="merged":
        df=df[df["cell_line"]==cell_type].reset_index(drop=True)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{cell_type}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{cell_type}.csv",index=False)



def write_ci_by_re(df_orig,re_type):
    df=df_orig.copy()
    df=df[df["re"]==re_type].reset_index(drop=True)
    df=calculate_tf_pair_cooperativity_index(df)
    df.to_csv(f"tf_pair_cooperativity_index_{re_type}.csv",index=False)
    df_tf=calculate_tf_cooperativity_index(df)
    df_tf.to_csv(f"tf_cooperativity_index_{re_type}.csv",index=False)



def write_tfs_codependent_and_redundant(path,suffix):
    df_tf=pd.read_csv(path)
    df_tf=df_tf[df_tf["c_sum"]>5].reset_index(drop=True)
    # redundant_tfs have cooperativity_index<0.3
    redundant_tfs=df_tf[df_tf["cooperativity_index"]<0.3]["protein2"].to_list()
    # sort alphabetically
    redundant_tfs.sort()
    with open(f"tfs_redundant_{suffix}.txt","w") as f:
        f.write("\n".join(redundant_tfs))
    codependent_tfs=df_tf[df_tf["cooperativity_index"]>0.7]["protein2"].to_list()
    codependent_tfs.sort()
    with open(f"tfs_codependent_{suffix}.txt","w") as f:
        f.write("\n".join(codependent_tfs))




df=read_all_files(0.1)
df["cell_line"]=df["dataset"].apply(lambda x: x.split("_")[-1])
for suffix in ["hepg2","k562","merged"]:
    write_ci_by_cell_type(df,suffix)




df["re"]=["promoter" if "promoter" in x else "enhancer" for x in df["dataset"]]
for suffix in ["promoter","enhancer"]:
    write_ci_by_re(df,suffix)


write_tfs_codependent_and_redundant("tf_cooperativity_index_merged.csv","merged")
write_tfs_codependent_and_redundant("tf_cooperativity_index_hepg2.csv","hepg2")
write_tfs_codependent_and_redundant("tf_cooperativity_index_k562.csv","k562")








# # ----------------------------------------------------
# # 4. The more codependent, the more picky about its parner
# # ----------------------------------------------------

# group by protein2 calculate entropy of cooperativity_index
from scipy.stats import entropy
df=pd.read_csv("tf_pair_cooperativity_index_post_filter.csv")

# group by protein2, calculate entropy of cooperativity_index, and count size of each group
df=df.groupby("protein2").agg({"cooperativity_index":entropy,"protein1":"count"}).reset_index()
# how many TFs 

df.rename(columns={"cooperativity_index":"entropy","protein1":"partner_count"},inplace=True)

# TF cooperativity ratio anti-correlates with entropy of cooperativity ratio
df_tf=pd.read_csv("tf_cooperativity_index.csv")
df=pd.merge(df,df_tf[["protein2","cooperativity_index"]],left_on="protein2",right_on="protein2",how="left")
df.dropna(inplace=True)
pearsonr(df["entropy"],df["cooperativity_index"]) 
spearmanr(df["entropy"],df["cooperativity_index"]) 
pearsonr(df["partner_count"],df["cooperativity_index"]) 
spearmanr(df["partner_count"],df["cooperativity_index"])





# ----------------------------------------------------
# 5. Redundant TFs have more versatile behavior in tf pairs
# ----------------------------------------------------
df_tf=pd.read_csv("tf_cooperativity_index.csv")
df_tf.rename(columns={"cooperativity_index":"cr_tf"},inplace=True)
df=pd.read_csv("tf_pair_cooperativity_index_post_filter.csv")
# df: group by protein2, calculate std cooperativity_index
df_tf_with_std=df.groupby("protein2").agg({"cooperativity_index":"std"}).reset_index()
df_tf_with_std.fillna(0,inplace=True)
df_tf_with_std.rename(columns={"cooperativity_index":"std"},inplace=True)
df_tf=df_tf.merge(df_tf_with_std,on="protein2",how="inner")
pearsonr(df_tf["cr_tf"],df_tf["std"]) # (-0.40,0)

# merge df_tf with df by protein2
df=pd.merge(df,df_tf[["protein2","cr_tf"]],on="protein2",how="left")
df=bin_and_label(df,"cr_tf",[0,0.2,0.4,0.6,0.8,1])
# violin plot
sns.violinplot(x="Bin",y="cooperativity_index",data=df)
plt.title("Distribution of TF pair cooperativity index")
plt.xlabel("TF cooperativity index")
plt.ylabel("TF pair cooperativity index")
# rotate x-axis labels by 45 degrees
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.3)
plt.savefig("Plots/tf_pair_cooperativity_index_by_tf_cooperativity_index.pdf")
plt.close()




# ----------------------------------------------------
# 3. Do TF pair relationship change across RE? 
# 85% TF pair show consistent cooperativity across sequence and cellular contexts
# but about 15% change from redundancy (promoter) to codependency(enhancer)
# ----------------------------------------------------
df=read_all_files()

df_promoters_hepg2=calculate_tf_pair_cooperativity_index(df[df["dataset"]=="promoter_hepg2"].reset_index(drop=True),"_promoters_hepg2")
df_promoters_k562=calculate_tf_pair_cooperativity_index(df[df["dataset"]=="promoter_k562"].reset_index(drop=True),"_promoters_k562")
df_enhancers_hepg2=calculate_tf_pair_cooperativity_index(df[df["dataset"]=="enhancer_hepg2"].reset_index(drop=True),"_enhancers_hepg2")
df_enhancers_k562=calculate_tf_pair_cooperativity_index(df[df["dataset"]=="enhancer_k562"].reset_index(drop=True),"_enhancers_k562")

# retain rows with sum_ci>0.5
df_promoters_hepg2=df_promoters_hepg2[df_promoters_hepg2["sum_ci_promoters_hepg2"]>0.5].reset_index(drop=True)
df_promoters_k562=df_promoters_k562[df_promoters_k562["sum_ci_promoters_k562"]>0.5].reset_index(drop=True)
df_enhancers_hepg2=df_enhancers_hepg2[df_enhancers_hepg2["sum_ci_enhancers_hepg2"]>0.5].reset_index(drop=True)
df_enhancers_k562=df_enhancers_k562[df_enhancers_k562["sum_ci_enhancers_k562"]>0.5].reset_index(drop=True)


# merge by protein1 and protein2, get all tf pairs from all 4 datasets
df_tf_pair=pd.merge(df_promoters_hepg2,df_promoters_k562,on=["protein1","protein2"],how="outer")
df_tf_pair=pd.merge(df_tf_pair,df_enhancers_hepg2,on=["protein1","protein2"],how="outer")
df_tf_pair=pd.merge(df_tf_pair,df_enhancers_k562,on=["protein1","protein2"],how="outer")
# for all column names with cooperativity_index, use "cr" instead
df_tf_pair.columns=[col.replace("cooperativity_index","cr") for col in df_tf_pair.columns]


# count redundancy case per row: cr<0.5
df_tf_pair["redundancy_case"]=0
df_tf_pair.loc[(df_tf_pair["cr_promoters_hepg2"]<0.5),"redundancy_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_promoters_k562"]<0.5),"redundancy_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_enhancers_hepg2"]<0.5),"redundancy_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_enhancers_k562"]<0.5),"redundancy_case"]+=1

# count codependency case per row: cr>0.5
df_tf_pair["codependency_case"]=0
df_tf_pair.loc[(df_tf_pair["cr_promoters_hepg2"]>0.5),"codependency_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_promoters_k562"]>0.5),"codependency_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_enhancers_hepg2"]>0.5),"codependency_case"]+=1
df_tf_pair.loc[(df_tf_pair["cr_enhancers_k562"]>0.5),"codependency_case"]+=1


# get confusing cases
df_confusing=df_tf_pair[(df_tf_pair["redundancy_case"]>0) & (df_tf_pair["codependency_case"]>0)].copy().reset_index(drop=True) # 1374/9140
# fill nan with 0
df_confusing.fillna(0,inplace=True)
# Does difference between enhancers and promoters matter cause majority of confusing cases?
# calculate cr for promoters
df_confusing["codependency_ci_promoters"]=df_confusing["codependency_ci_promoters_hepg2"]+df_confusing["codependency_ci_promoters_k562"]
df_confusing["redundancy_ci_promoters"]=df_confusing["redundancy_ci_promoters_hepg2"]+df_confusing["redundancy_ci_promoters_k562"]
df_confusing["sum_ci_promoters"]=df_confusing["redundancy_ci_promoters"].abs()+df_confusing["codependency_ci_promoters"]
df_confusing["cr_promoters"]=df_confusing["codependency_ci_promoters"]/(df_confusing["redundancy_ci_promoters"].abs()+df_confusing["codependency_ci_promoters"])
# remove columns
df_confusing.drop(["codependency_ci_promoters_hepg2","redundancy_ci_promoters_hepg2","cr_promoters_hepg2","cr_promoters_k562",'sum_ci_promoters_hepg2','sum_ci_promoters_k562',
                   "codependency_ci_promoters_k562","redundancy_ci_promoters_k562"],axis=1,inplace=True)

# calculate cr for enhancers
df_confusing["codependency_ci_enhancers"]=df_confusing["codependency_ci_enhancers_hepg2"]+df_confusing["codependency_ci_enhancers_k562"]
df_confusing["redundancy_ci_enhancers"]=df_confusing["redundancy_ci_enhancers_hepg2"]+df_confusing["redundancy_ci_enhancers_k562"]
df_confusing["sum_ci_enhancers"]=df_confusing["redundancy_ci_enhancers"].abs()+df_confusing["codependency_ci_enhancers"]
df_confusing["cr_enhancers"]=df_confusing["codependency_ci_enhancers"]/(df_confusing["redundancy_ci_enhancers"].abs()+df_confusing["codependency_ci_enhancers"])
# remove columns
df_confusing.drop(["codependency_ci_enhancers_hepg2","redundancy_ci_enhancers_hepg2","cr_enhancers_hepg2","cr_enhancers_k562",'sum_ci_enhancers_hepg2','sum_ci_enhancers_k562',
                   "codependency_ci_enhancers_k562","redundancy_ci_enhancers_k562"],axis=1,inplace=True)

df_confusing[["cr_promoters","cr_enhancers"]]
# count number of rows with Nan
df_confusing["cr_promoters"].isna().sum() # 81
# count number of rows where cr_promoters<0.5 and cr_enhancers>0.5
df_confusing[(df_confusing["cr_promoters"]<0.5) & (df_confusing["cr_enhancers"]>0.5)] # 3727
df_confusing[(df_confusing["cr_promoters"]>0.5) & (df_confusing["cr_enhancers"]<0.5)] # 291







# ----------------------------------------------------
# 4. Do TF cooperativity change across RE? 
# ----------------------------------------------------


df_tf[df_tf["cr_promoters"]<df_tf["cr_enhancers"]] # 318

# which tf have "cr_promoters"<0.3 and "cr_enhancers">0.7
df_tf[(df_tf["cr_promoters"]<0.4) & (df_tf["cr_enhancers"]>0.7)]

df_tf[df_tf["cr_enhancers"]==1]