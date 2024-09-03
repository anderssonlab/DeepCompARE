import pandas as pd
from loguru import logger
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr,spearmanr
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import read_cooperativity


# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def analysis(file_prefix,sep):
    if "hepg2" in file_prefix:
        track_nums=[0,2,4,6]
    if "k562" in file_prefix:
        track_nums=[1,3,5,7]
    if "common" in file_prefix:
        track_nums=[0,1,2,3,4,5,6,7]
    df=read_cooperativity(f"/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_{file_prefix}.csv",track_nums=track_nums)
    # group by region_idx and codependency, count the number of zero and one in codependency column in each group
    df_grouped=df.groupby(["region_idx","codependency"])[["codependency"]].count().unstack().fillna(0).reset_index()
    # reduce to single index
    df_grouped.columns = df_grouped.columns.droplevel(0)
    df_grouped.columns=["region_idx","redundant","codependent"]
    df_confusing_re=df_grouped[(df_grouped["redundant"]>0) & (df_grouped["codependent"]>0)]
    logger.info(f"# Total RE regions: {df_grouped.shape[0]}")
    logger.info(f"# Confusing RE regions (contradiction within RE): {df_confusing_re.shape[0]}")
    # add region information
    regions_df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{file_prefix}.bed",sep=sep,header=None)
    # label promoters i with f"Seq{i}"
    regions_df["region_idx"]=regions_df.index    
    regions_df["region_idx"]=regions_df["region_idx"].apply(lambda x: f"Region{x}")
    # merge with df_grouped
    df_grouped=pd.merge(df_grouped,regions_df,on="region_idx",how="left")
    df_grouped.to_csv(f"{file_prefix}_tf_pair_cooperativity_counts_lenient.csv",index=False)


# ----------------------------------------------------
# 1. get promoter/enhancer level TF pair info
# ----------------------------------------------------
analysis("enhancers_hepg2","\t")
analysis("enhancers_k562","\t")
analysis("promoters_hepg2","\t")
analysis("promoters_k562","\t")


df_enhancers_hepg2=pd.read_csv("enhancers_hepg2_tf_pair_cooperativity_counts_lenient.csv")
df_enhancers_k562=pd.read_csv("enhancers_k562_tf_pair_cooperativity_counts_lenient.csv")
df_promoters_hepg2=pd.read_csv("promoters_hepg2_tf_pair_cooperativity_counts_lenient.csv")
df_promoters_k562=pd.read_csv("promoters_k562_tf_pair_cooperativity_counts_lenient.csv")

# count rows with codependent=0 or redundant=0
df_sub=df_promoters_k562
df_sub.shape
df_sub[(df_sub["codependent"]==0)].shape[0]
df_sub[(df_sub["redundant"]==0)].shape[0]



df_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_enhancers_k562["dataset"]="enhancers_k562"
df_promoters_hepg2["dataset"]="promoters_hepg2"
df_promoters_k562["dataset"]="promoters_k562"

df_enhancers_hepg2["codependent_pair_percentage"]=df_enhancers_hepg2["codependent"]/(df_enhancers_hepg2["codependent"]+df_enhancers_hepg2["redundant"])
df_enhancers_k562["codependent_pair_percentage"]=df_enhancers_k562["codependent"]/(df_enhancers_k562["codependent"]+df_enhancers_k562["redundant"])
df_promoters_hepg2["codependent_pair_percentage"]=df_promoters_hepg2["codependent"]/(df_promoters_hepg2["codependent"]+df_promoters_hepg2["redundant"])
df_promoters_k562["codependent_pair_percentage"]=df_promoters_k562["codependent"]/(df_promoters_k562["codependent"]+df_promoters_k562["redundant"])


df=pd.concat([df_enhancers_hepg2,df_enhancers_k562,df_promoters_hepg2,df_promoters_k562],axis=0)


stat,p_hepg2=mannwhitneyu(df_enhancers_hepg2["codependent_pair_percentage"],df_promoters_hepg2["codependent_pair_percentage"])
df_enhancers_hepg2["codependent_pair_percentage"].median()
df_promoters_hepg2["codependent_pair_percentage"].median()

stat,p_k562=mannwhitneyu(df_enhancers_k562["codependent_pair_percentage"],df_promoters_k562["codependent_pair_percentage"])
df_enhancers_k562["codependent_pair_percentage"].median()
df_promoters_k562["codependent_pair_percentage"].median()



# kdeplot codependent_pair_percentage, hue=dataset
sns.kdeplot(data=df,x="codependent_pair_percentage",hue="dataset",common_norm=False)
plt.title("Distribution of codependent pair percentage")
plt.savefig("codependent_pair_percentage_lenient.png")
plt.close()




# ----------------------------------------------------
# 2. Compare between sharp and broad promoters
# ----------------------------------------------------


analysis("promoters_broad_common"," ")
analysis("promoters_sharp_common"," ")


analysis("promoters_broad_hepg2"," ")
analysis("promoters_sharp_hepg2"," ")
analysis("promoters_broad_k562"," ")
analysis("promoters_sharp_k562"," ")


# calculate sub ratio distribution for sharp and broad
df_broad1=pd.read_csv("promoters_broad_common_tf_pair_cooperativity_counts_lenient.csv")
df_broad2=pd.read_csv("promoters_broad_hepg2_tf_pair_cooperativity_counts_lenient.csv")
df_broad3=pd.read_csv("promoters_broad_k562_tf_pair_cooperativity_counts_lenient.csv")
df_sharp1=pd.read_csv("promoters_sharp_common_tf_pair_cooperativity_counts_lenient.csv")
df_sharp2=pd.read_csv("promoters_sharp_hepg2_tf_pair_cooperativity_counts_lenient.csv")
df_sharp3=pd.read_csv("promoters_sharp_k562_tf_pair_cooperativity_counts_lenient.csv")

df_broad1["dataset"]="broad_common"
df_broad2["dataset"]="broad_hepg2"
df_broad3["dataset"]="broad_k562"
df_sharp1["dataset"]="sharp_common"
df_sharp2["dataset"]="sharp_hepg2"
df_sharp3["dataset"]="sharp_k562"

df=pd.concat([df_broad1,df_broad2,df_broad3,df_sharp1,df_sharp2,df_sharp3],axis=0)
df["promoter_type"]=df["dataset"].apply(lambda x: x.split("_")[0])
df["dataset"]=df["dataset"].apply(lambda x: x.split("_")[1])


df["total"]=df["redundant"]+df["codependent"]
df["codepedent_pair_percentage"]=df["codependent"]/df["total"]


mannwhitneyu(df[df["promoter_type"]=="broad"]["codepedent_pair_percentage"],df[df["promoter_type"]=="sharp"]["codepedent_pair_percentage"])

df[df["promoter_type"]=="broad"]["codepedent_pair_percentage"].median()
df[df["promoter_type"]=="sharp"]["codepedent_pair_percentage"].median()

# violinplot codependent_pair_percentage, hue=dataset
sns.kdeplot(data=df,x="codepedent_pair_percentage",hue="promoter_type",common_norm=False)
# add annotation, p round to 2 decimal
plt.title("Distribution of codependent pair percentage")
plt.savefig("codependent_pair_percentage_sharp_broad.png")
plt.close()




# calculate pearson correlation and spearman correlation 
pearsonr(df["codepedent_pair_percentage"],df["14"])
spearmanr(df["codepedent_pair_percentage"],df["14"])
# reverse order???

# ----------------------------------------------------
# 3. Why do promoter and enhancer differ in additivity profile?
# H1: same TF pair has different additivity profile 
# H2: different TF pair distribution
# ----------------------------------------------------

cell_type="hepg2"
remap="Hep-G2"

cell_type="k562"
remap="K-562"

df_promoter=aggregate_tracks("promoters",cell_type,remap)
# group by protein1, protein2, sum sub_additive and super_additive
df_promoter=df_promoter.groupby(["protein1","protein2"])[["sub_additive","super_additive"]].sum().reset_index()
df_promoter.rename(columns={"sub_additive":"sub_additive_promoter","super_additive":"super_additive_promoter"},inplace=True)
df_enhancer=aggregate_tracks("enhancers",cell_type,remap)
df_enhancer=df_enhancer.groupby(["protein1","protein2"])[["sub_additive","super_additive"]].sum().reset_index()
df_enhancer.rename(columns={"sub_additive":"sub_additive_enhancer","super_additive":"super_additive_enhancer"},inplace=True)

# merge promoter and enhancer by 'region_idx', 'protein1', 'protein2', 'chromosome1', 'start_rel1','end_rel1', 'strand', 'score', 'chromosome2', 'start_rel2', 'end_rel2', 'score2'
df=pd.merge(df_promoter,df_enhancer,on=['protein1','protein2'],how="outer")
# replace NaN with 0
df.fillna(0,inplace=True)


df["additivity_promoter"]="unknown"
df.loc[df["sub_additive_promoter"]>df["super_additive_promoter"],"additivity_promoter"]="sub"
df.loc[df["sub_additive_promoter"]<df["super_additive_promoter"],"additivity_promoter"]="super"

df["additivity_enhancer"]="unknown"
df.loc[df["sub_additive_enhancer"]>df["super_additive_enhancer"],"additivity_enhancer"]="sub"
df.loc[df["sub_additive_enhancer"]<df["super_additive_enhancer"],"additivity_enhancer"]="super"

df_known=df[(df["additivity_promoter"]!="unknown") & (df["additivity_enhancer"]!="unknown")].reset_index(drop=True)
# how many TF pairs have different additivity profile
df_diff=df_known[df_known["additivity_promoter"]!=df_known["additivity_enhancer"]].reset_index(drop=True)



# Conclusions: 
# For HepG2 
# union(tf_pair_promoter,tf_pair_enhancer)=7068
# intersection(tf_pair_promoter,tf_pair_enhancer)=2264
# inconsistent TF_pair behavior between promoter and enhancer: 802

# For K562:
# union(tf_pair_promoter,tf_pair_enhancer)=10804
# intersection(tf_pair_promoter,tf_pair_enhancer)=4050
# inconsistent TF_pair behavior between promoter and enhancer: 1320

