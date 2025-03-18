import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

from scipy.stats import pearsonr



re="distal"


# Major confounder: Only investigated high quality TFBS, not bound TFBSs
# TODO: subset using chip seq to contain only bound regions


def count_mediator_interactors_per_region(df_med,mediator,df_tfbs):
    if mediator is not "all":
        df_med=df_med[df_med["bait"]==mediator].reset_index(drop=True)
    # is protein in df_med.gene
    df=df_tfbs.copy()
    df["mediator"]=df["protein"].apply(lambda x: x in df_med["gene"].values.tolist())
    df["mediator"]=df["mediator"].astype(int)
    # group by region, and count the number of mediator, get the mean gc content
    df=df.groupby("region").agg({"mediator":"sum"}).reset_index()
    return df



def analyze_tfbs(df_histone, df_tfbs, df_med, mediator, re):
    # count number of mediator-interactor TFBSs per region
    df_res = count_mediator_interactors_per_region(df_med, mediator, df_tfbs)
    df_res = pd.merge(df_res, df_histone, on="region", how="inner")
    # plot distribution of columns starting with "log_signal"
    log_signal_cols = [col for col in df_histone.columns if "log_signal" in col]
    log_signal_cols = [col.replace("log_signal_", "") for col in log_signal_cols]
    for col in log_signal_cols:
        logger.info(f"plotting {col}")
        plt.figure(figsize=(3, 2.5))
        # calculate pearson correlation and p-value
        r, p = pearsonr(df_res["mediator"], df_res["log_signal_" + col])
        sns.boxplot(x="mediator", y="log_signal_" + col, data=df_res, fliersize=0, linewidth=0.5)
        # Count occurrences of each x category
        category_counts = df_res["mediator"].value_counts().sort_index()
        # Modify x-tick labels to include counts
        xtick_labels = [f"{x} (n={category_counts[x]})" for x in category_counts.index]
        plt.xticks(ticks=range(len(category_counts)), labels=xtick_labels,rotation=45, fontsize=5)
        plt.yticks(fontsize=5)
        plt.title(f"{col}, {re}, r={r:.2f}, p={p:.2f}", fontsize=7)
        plt.xlabel(f"# {mediator}-interactors TFBSs",fontsize=7)
        plt.ylabel(f"log {col}",fontsize=7)
        plt.tight_layout()
        plt.savefig(f"{col}_{mediator}_{re}.pdf")
        plt.close()
        
        
        
        



# read histone marks
df_histone=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/{re}_k562.csv")
# TODO: choose to -1 or not
# df_histone["end"]=df_histone["end"]-1
# add region in format chr:start-end
df_histone["region"]=df_histone["chrom"].astype(str)+":"+df_histone["start"].astype(str)+"-"+df_histone["end"].astype(str)



# plot correlation between histone marks
corr_map=df_histone.iloc[:,3:].corr()
# rename columns: remove log_signal_
corr_map.columns=[col.replace("log_signal_","") for col in corr_map.columns]
corr_map.index=[col.replace("log_signal_","") for col in corr_map.index]
# plot clustermap
plt.figure()
sns.clustermap(corr_map)
plt.savefig(f"clustermap_histone_{re}.pdf")
plt.close()




# read tfbs
df_tfbs=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_dhs_{re}_k562.csv")
# subset for proper background TFs
proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
df_tfbs=df_tfbs[df_tfbs['protein'].isin(proteins_background)].reset_index(drop=True)


# read mediator interactors
df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)


# analyze # mediator interactors v.s. histone marks
for mediator in df_med["bait"].unique().tolist()+["all"]:
    analyze_tfbs(df_histone,df_tfbs,df_med,mediator,re)
    logger.info(f"{mediator} done")
    
    










#------------------------------------------------
#  avg number of TFBSs per promoter v.s. enhancer
#------------------------------------------------

# df_promoter=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_k562.csv")
# df_enhancer=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_k562.csv")

# df_med=pd.read_csv(f"2024-08-07_MED-TF_interactions.txt",sep="\t")
# df_med=df_med[df_med["significant"]].reset_index(drop=True)
# tfs=df_med["gene"].unique().tolist()

# # select only tfs
# df_promoter=df_promoter[df_promoter["protein"].isin(tfs)].reset_index(drop=True)
# df_enhancer=df_enhancer[df_enhancer["protein"].isin(tfs)].reset_index(drop=True)

# tf_counts_promoter=df_promoter.protein.value_counts()
# tf_counts_enhancer=df_enhancer.protein.value_counts()

# # how many promoter regions in total?
# n_promoter_regions=len(df_promoter["region"].unique())
# n_enhancer_regions=len(df_enhancer["region"].unique())

# # normalize by number of regions
# tf_counts_promoter=tf_counts_promoter/n_promoter_regions
# tf_counts_enhancer=tf_counts_enhancer/n_enhancer_regions

# # merge tf counts by name 
# df=pd.merge(tf_counts_promoter,tf_counts_enhancer,left_index=True,right_index=True,how="inner",suffixes=("_promoter","_enhancer"))

# # rename protein_promoter to promoter
# df.rename(columns={"protein_promoter":"promoter","protein_enhancer":"enhancer"},inplace=True)

# # scatter plot 
# plt.figure()
# plt.scatter(df["promoter"],df["enhancer"])
# # add name for each point
# for i in range(len(df)):
#     plt.text(df["promoter"][i],df["enhancer"][i],df.index[i])

# # add diagonal line
# min_val=min(df["promoter"].min(),df["enhancer"].min())
# max_val=max(df["promoter"].max(),df["enhancer"].max())
# plt.plot([min_val,max_val],[min_val,max_val],color="gray",linestyle="--",linewidth=0.5)
# plt.xlabel("avg # TFBSs per promoter region")
# plt.ylabel("avg # TFBSs per enhancer region")
# plt.title("avg # TFBSs per region")
# plt.savefig("avg_tfbs_per_region.png")
# plt.close()














# nohup python3 tfbs.py > tfbs.out &