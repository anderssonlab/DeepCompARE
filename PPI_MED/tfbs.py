import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

re="promoters"


# Major confounder: Only investigated high quality TFBS, not bound TFBSs
# TODO: subset using chip seq to contain only bound regions
# TODO: annotate number. group minority groups
# TODO: add pearson r and p
# TODO: overall correlation between the histone modifications, is corr(ATAC,H2AC)<0

def count_mediator_interactors_per_region(df_med,mediator,df_tfbs):
    if mediator is not "all":
        df_med=df_med[df_med["bait"]==mediator].reset_index(drop=True)
    # is protein in df_med.gene
    df=df_tfbs.copy()
    df["mediator"]=df["protein"].apply(lambda x: x in df_med["gene"].values.tolist())
    df["mediator"]=df["mediator"].astype(int)
    # group by region, and count the number of mediator
    df=df.groupby("region").agg({"mediator":"sum"}).reset_index()
    return df



def analyze_tfbs(df_histone,df_tfbs,df_med,mediator,re):
    # count number of mediator-interactor TFBSs per region
    df_res=count_mediator_interactors_per_region(df_med,mediator,df_tfbs)
    df_res=pd.merge(df_res,df_histone,on="region",how="inner")
    # plot distribution of columns starting with "log_signal"
    log_signal_cols=[col for col in df_histone.columns if "log_signal" in col]
    # remove "log_signal" from column names, for shorter output file names
    log_signal_cols=[col.replace("log_signal_","") for col in log_signal_cols]
    # plot all histone marks
    for col in log_signal_cols:
        logger.info(f"plotting {col}")
        plt.figure()
        sns.boxplot(x="mediator",y="log_signal_"+col,data=df_res,fliersize=0)
        plt.title(f"{col} signal in {re}")
        plt.xlabel(f"# {mediator}-interactor TFBSs")
        plt.savefig(f"{col}_{mediator}_{re}.png")
        plt.close()


df_histone=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/{re}_k562.csv")
df_histone["end"]=df_histone["end"]-1
# add region in format chr:start-end
df_histone["region"]=df_histone["chrom"].astype(str)+":"+df_histone["start"].astype(str)+"-"+df_histone["end"].astype(str)

df_tfbs=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{re}_k562.csv")
# subset for proper background TFs
proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
df_tfbs=df_tfbs[df_tfbs['protein'].isin(proteins_background)].reset_index(drop=True)

df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)
for mediator in df_med["bait"].unique().tolist()+["all"]:
    analyze_tfbs(df_histone,df_tfbs,df_med,mediator,re)
    logger.info(f"{mediator} done")
    
    











#------------------------------------------------
#  avg number of TFBSs per promoter v.s. enhancer
#------------------------------------------------

df_promoter=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_promoters_k562.csv")
df_enhancer=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_enhancers_k562.csv")

df_med=pd.read_csv(f"2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)
tfs=df_med["gene"].unique().tolist()

# select only tfs
df_promoter=df_promoter[df_promoter["protein"].isin(tfs)].reset_index(drop=True)
df_enhancer=df_enhancer[df_enhancer["protein"].isin(tfs)].reset_index(drop=True)

tf_counts_promoter=df_promoter.protein.value_counts()
tf_counts_enhancer=df_enhancer.protein.value_counts()

# how many promoter regions in total?
n_promoter_regions=len(df_promoter["region"].unique())
n_enhancer_regions=len(df_enhancer["region"].unique())

# normalize by number of regions
tf_counts_promoter=tf_counts_promoter/n_promoter_regions
tf_counts_enhancer=tf_counts_enhancer/n_enhancer_regions

# merge tf counts by name 
df=pd.merge(tf_counts_promoter,tf_counts_enhancer,left_index=True,right_index=True,how="inner",suffixes=("_promoter","_enhancer"))

# rename protein_promoter to promoter
df.rename(columns={"protein_promoter":"promoter","protein_enhancer":"enhancer"},inplace=True)

# scatter plot 
plt.figure()
plt.scatter(df["promoter"],df["enhancer"])
# add name for each point
for i in range(len(df)):
    plt.text(df["promoter"][i],df["enhancer"][i],df.index[i])

# add diagonal line
min_val=min(df["promoter"].min(),df["enhancer"].min())
max_val=max(df["promoter"].max(),df["enhancer"].max())
plt.plot([min_val,max_val],[min_val,max_val],color="gray",linestyle="--",linewidth=0.5)
plt.xlabel("avg # TFBSs per promoter region")
plt.ylabel("avg # TFBSs per enhancer region")
plt.title("avg # TFBSs per region")
plt.savefig("avg_tfbs_per_region.png")
plt.close()














# nohup python3 tfbs.py > tfbs.out &