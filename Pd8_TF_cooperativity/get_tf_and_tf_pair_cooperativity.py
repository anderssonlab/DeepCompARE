import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,spearmanr
from adjustText import adjust_text
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from utils import split_dimer, bin_and_label
from tf_cooperativity import read_cooperativity, calculate_tf_pair_cooperativity_ratio



# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def read_all_files():
    # read cooperativity
    # TODO: choose whether to use lenient or strict
    df_promoter_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_hepg2.csv",track_nums=[0,2,4,6])
    df_promoter_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_promoters_k562.csv",track_nums=[1,3,5,7])
    df_enhancer_hepg2=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_hepg2.csv",track_nums=[0,2,4,6])
    df_enhancer_k562=read_cooperativity("/isdata/alab/people/pcr980/DeepCompare/Pd6_mutate_pair/mutate_pairs_lenient_enhancers_k562.csv",track_nums=[1,3,5,7])
    df_promoter_hepg2["dataset"]="promoter_hepg2"
    df_promoter_k562["dataset"]="promoter_k562"
    df_enhancer_hepg2["dataset"]="enhancer_hepg2"
    df_enhancer_k562["dataset"]="enhancer_k562"
    df=pd.concat([df_promoter_hepg2,df_promoter_k562,df_enhancer_hepg2,df_enhancer_k562],axis=0)
    return df



# ----------------------------------------------------
# 1. Get TF pair cooperativity ratio
# ----------------------------------------------------
df=read_all_files()
df=calculate_tf_pair_cooperativity_ratio(df)
df.to_csv("tf_pair_cooperativity_ratio_pre_filter_lenient.csv",index=False)

# select sum_ci>1
df["sum_ci"].describe()
df=df[df["sum_ci"]>1].reset_index(drop=True)
df.to_csv("tf_pair_cooperativity_ratio_post_filter_lenient.csv",index=False)


# ----------------------------------------------------
# 2. Plot TF pair cooperativity ratio and show extreme pairs
# ----------------------------------------------------
# use post_filter for plotting
df=pd.read_csv("tf_pair_cooperativity_ratio_post_filter_lenient.csv")
# how many TF pairs have cooperativity_ratio between 0.3 and 0.7
df[(df["cooperativity_ratio"]>0.3) & (df["cooperativity_ratio"]<0.7)].shape[0] 


# pearson correlation
pearsonr(df["redundancy_ci"],df["codependency_ci"]) # (0.55,0)
pearsonr(df["sum_ci"],df["cooperativity_ratio"]) 

# plot histogram
sns.histplot(df["cooperativity_ratio"],kde=True)
plt.title("TF pair cooperativity atio distribution")
plt.savefig("Plots/tf_pair_cooperativity_ratio_distribution_lenient.png")
plt.close()

# plot heatmap
df=df.pivot(index="protein1",columns="protein2",values="cooperativity_ratio")
sns.heatmap(df,cmap="coolwarm",vmin=0,vmax=1)
plt.title("TF pair cooperativity ratio")
plt.subplots_adjust(bottom=0.3)
plt.subplots_adjust(left=0.3)
plt.savefig("Plots/tf_pair_cooperativity_ratio_heatmap_lenient.png")
plt.close()# plot heatmap


# ----------------------------------------------------
# 3. Get TF cooperativity ratio and histogran
# ----------------------------------------------------
# use pre_filter for TF level aggregation
df=pd.read_csv("tf_pair_cooperativity_ratio_pre_filter_lenient.csv")
df_tf=df.groupby("protein2").agg({"redundancy_ci":"sum",
                                  "codependency_ci":"sum",
                                  "cooperativity_ratio":"std"}).reset_index()
df_tf["sum_ci"]=df_tf["redundancy_ci"].abs()+df_tf["codependency_ci"]
df_tf=df_tf[df_tf["sum_ci"]>5].reset_index(drop=True)
df_tf.rename(columns={"cooperativity_ratio":"std"},inplace=True)
df_tf["cooperativity_ratio"]=df_tf["codependency_ci"]/(df_tf["redundancy_ci"].abs()+df_tf["codependency_ci"])
df_tf.sort_values("cooperativity_ratio",ascending=False,inplace=True)
df_tf.to_csv("tf_cooperativity_ratio_lenient.csv",index=False)

sns.histplot(df_tf["cooperativity_ratio"],kde=True)
plt.title("TF cooperativity ratio distribution")
plt.savefig("Plots/tf_cooperativity_ratio_distribution_lenient.png")
plt.close()


# ----------------------------------------------------
# 4. The more codependent, the more picky about its parner
# ----------------------------------------------------

# group by protein2 calculate entropy of cooperativity_ratio
from scipy.stats import entropy
df=pd.read_csv("tf_pair_cooperativity_ratio_post_filter_lenient.csv")
# group by protein2, calculate entropy of cooperativity_ratio, and count size of each group
df=df.groupby("protein2").agg({"cooperativity_ratio":entropy,"protein1":"count"}).reset_index()

df.rename(columns={"cooperativity_ratio":"entropy","protein1":"partner_count"},inplace=True)

# TF cooperativity ratio anti-correlates with entropy of cooperativity ratio
df_tf=pd.read_csv("tf_cooperativity_ratio.csv")
df=pd.merge(df,df_tf[["protein2","cooperativity_ratio"]],left_on="protein2",right_on="protein2",how="left")
df.dropna(inplace=True)
pearsonr(df["entropy"],df["cooperativity_ratio"]) 
spearmanr(df["entropy"],df["cooperativity_ratio"]) 
pearsonr(df["partner_count"],df["cooperativity_ratio"]) 
spearmanr(df["partner_count"],df["cooperativity_ratio"])

# ----------------------------------------------------
# 4. Get tfs_redundant and tfs_codependent
# ----------------------------------------------------
df_tf=pd.read_csv("tf_cooperativity_ratio_lenient.csv")

# redundant_tfs have cooperativity_ratio<0.25
redundant_tfs=df_tf[df_tf["cooperativity_ratio"]<0.25]["protein2"].to_list()
redundant_tfs=[tf for tf in redundant_tfs if "::" not in tf]
with open("tfs_redundant_lenient.txt","w") as f:
    f.write("\n".join(redundant_tfs))


codependent_tfs=df_tf[df_tf["cooperativity_ratio"]>0.75]["protein2"].to_list()
with open("tfs_codependent_lenient.txt","w") as f:
    f.write("\n".join(split_dimer(codependent_tfs)))




# ----------------------------------------------------
# 5. Redundant TFs have more versatile behavior in tf pairs
# ----------------------------------------------------
df_tf=pd.read_csv("tf_cooperativity_ratio_lenient.csv")
df_tf.rename(columns={"cooperativity_ratio":"cr_tf"},inplace=True)
df=pd.read_csv("tf_pair_cooperativity_ratio_post_filter_lenient.csv")
pearsonr(df_tf["cr_tf"],df_tf["std"]) # (-0.32,0)

# merge df_tf with df by protein2
df=pd.merge(df,df_tf[["protein2","cr_tf"]],on="protein2",how="left")
df=bin_and_label(df,"cr_tf",[0,0.2,0.4,0.6,0.8,1])
# violin plot
sns.violinplot(x="Bin",y="cooperativity_ratio",data=df)
plt.title("Distribution of TF pair cooperativity ratio")
plt.xlabel("TF cooperativity ratio")
plt.ylabel("TF pair cooperativity ratio")
# rotate x-axis labels by 45 degrees
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.3)
plt.savefig("Plots/tf_pair_cooperativity_ratio_by_tf_cooperativity_ratio_lenient.png")
plt.close()




# ----------------------------------------------------
# 3. Do TF pair relationship change across RE? 
# 85% TF pair show consistent cooperativity across sequence and cellular contexts
# but about 10% change from redundancy (promoter) to codependency(enhancer)
# ----------------------------------------------------
df=read_all_files()

df_promoters_hepg2=calculate_tf_pair_cooperativity_ratio(df[df["dataset"]=="promoter_hepg2"].reset_index(drop=True),"_promoters_hepg2")
df_promoters_k562=calculate_tf_pair_cooperativity_ratio(df[df["dataset"]=="promoter_k562"].reset_index(drop=True),"_promoters_k562")
df_enhancers_hepg2=calculate_tf_pair_cooperativity_ratio(df[df["dataset"]=="enhancer_hepg2"].reset_index(drop=True),"_enhancers_hepg2")
df_enhancers_k562=calculate_tf_pair_cooperativity_ratio(df[df["dataset"]=="enhancer_k562"].reset_index(drop=True),"_enhancers_k562")

# retain rows with sum_ci>0.5
df_promoters_hepg2=df_promoters_hepg2[df_promoters_hepg2["sum_ci_promoters_hepg2"]>0.5].reset_index(drop=True)
df_promoters_k562=df_promoters_k562[df_promoters_k562["sum_ci_promoters_k562"]>0.5].reset_index(drop=True)
df_enhancers_hepg2=df_enhancers_hepg2[df_enhancers_hepg2["sum_ci_enhancers_hepg2"]>0.5].reset_index(drop=True)
df_enhancers_k562=df_enhancers_k562[df_enhancers_k562["sum_ci_enhancers_k562"]>0.5].reset_index(drop=True)


# merge by protein1 and protein2
df_tf_pair=pd.merge(df_promoters_hepg2,df_promoters_k562,on=["protein1","protein2"],how="outer")
df_tf_pair=pd.merge(df_tf_pair,df_enhancers_hepg2,on=["protein1","protein2"],how="outer")
df_tf_pair=pd.merge(df_tf_pair,df_enhancers_k562,on=["protein1","protein2"],how="outer")
# for all column names with cooperativity_ratio, use "cr" instead
df_tf_pair.columns=[col.replace("cooperativity_ratio","cr") for col in df_tf_pair.columns]

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
df_confusing[(df_confusing["cr_promoters"]<0.5) & (df_confusing["cr_enhancers"]>0.5)] # 938
df_confusing[(df_confusing["cr_promoters"]>0.5) & (df_confusing["cr_enhancers"]<0.5)] # 78







# ----------------------------------------------------
# 4. Do TF cooperativity change across RE? 
# ----------------------------------------------------
df=read_all_files()
df["re"]=["promoter" if "promoter" in x else "enhancer" for x in df["dataset"]]

df_promoters=calculate_tf_pair_cooperativity_ratio(df[df["re"]=="promoter"].reset_index(drop=True))
df_enhancers=calculate_tf_pair_cooperativity_ratio(df[df["re"]=="enhancer"].reset_index(drop=True))

# group by protein2, sum redundancy_ci, codependency_ci, sum_ci

df_promoters=df_promoters.groupby("protein2").agg({"redundancy_ci":"sum","codependency_ci":"sum"}).reset_index()
df_promoters["sum_ci"]=df_promoters["redundancy_ci"].abs()+df_promoters["codependency_ci"]
df_promoters=df_promoters[df_promoters["sum_ci"]>2].reset_index(drop=True)
df_promoters["cooperativity_ratio"]=df_promoters["codependency_ci"].abs()/(df_promoters["redundancy_ci"].abs()+df_promoters["codependency_ci"])

df_enhancers=df_enhancers.groupby("protein2").agg({"redundancy_ci":"sum","codependency_ci":"sum"}).reset_index()
df_enhancers["sum_ci"]=df_enhancers["redundancy_ci"].abs()+df_enhancers["codependency_ci"]
df_enhancers=df_enhancers[df_enhancers["sum_ci"]>2].reset_index(drop=True)
df_enhancers["cooperativity_ratio"]=df_enhancers["codependency_ci"].abs()/(df_enhancers["redundancy_ci"].abs()+df_enhancers["codependency_ci"])

df_promoters=df_promoters[["protein2","cooperativity_ratio"]].rename(columns={"cooperativity_ratio":"cr_promoters"})
df_enhancers=df_enhancers[["protein2","cooperativity_ratio"]].rename(columns={"cooperativity_ratio":"cr_enhancers"})
df_tf=pd.merge(df_promoters,df_enhancers,on=["protein2"],how="outer")
df_tf.to_csv("tf_cooperativity_ratio_promoters_vs_enhancers_lenient.csv",index=False)

# select rows with at least one Nan
df_tf[df_tf.isna().any(axis=1)] # 0
# remove rows with Nan
df_tf.dropna(inplace=True)
# add tf_type
tfs_redundant=pd.read_csv("tfs_redundant_lenient.txt",header=None)[0].to_list()
tfs_codependent=pd.read_csv("tfs_codependent_lenient.txt",header=None)[0].to_list()
df_tf["tf_type"]="undetermined"
df_tf.loc[df_tf["protein2"].isin(tfs_redundant),"tf_type"]="redundant"
df_tf.loc[df_tf["protein2"].isin(tfs_codependent),"tf_type"]="codependent"

# sns scatter plot,small dot
# adjust text
plt.figsize=(8,9)
sns.scatterplot(x="cr_promoters",y="cr_enhancers",data=df_tf,hue="tf_type",s=5)
plt.title("TF cooperativity ratio in promoters and enhancers")
plt.xlabel("TF cooperativity ratio in promoters")
plt.ylabel("TF cooperativity ratio in enhancers")
# add abline
plt.plot([0,1],[0,1],color="black",linestyle="--")
# annotate offdiagonal points
# adjust text to avoid overlap,small text
texts=[]
for i in range(df_tf.shape[0]):
    if abs(df_tf.iloc[i]["cr_promoters"]-df_tf.iloc[i]["cr_enhancers"])>0.4:
        texts.append(plt.text(df_tf.iloc[i]["cr_promoters"],df_tf.iloc[i]["cr_enhancers"],df_tf.iloc[i]["protein2"],fontsize=6))
    if df_tf.iloc[i]["cr_promoters"]>df_tf.iloc[i]["cr_enhancers"]:
        texts.append(plt.text(df_tf.iloc[i]["cr_promoters"],df_tf.iloc[i]["cr_enhancers"],df_tf.iloc[i]["protein2"],fontsize=6))

adjust_text(texts,arrowprops=dict(arrowstyle="-",color='black'))
# make legend at bottom right
plt.legend(loc='lower right')
plt.savefig("Plots/tf_cooperativity_ratio_promoters_vs_enhancers_lenient.png")
plt.close()

