import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr

from loguru import logger


re="promoters"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{re}_k562.csv")
# subset for proper background TFs
proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
df=df[df['protein'].isin(proteins_background)].reset_index(drop=True)



df["diff_activity_probability_cage"]=df["isa_track1"]-df["isa_track9"]
df["diff_activity_probability_dhs"]=df["isa_track3"]-df["isa_track11"]
df["diff_activity_probability_starr"]=df["isa_track5"]-df["isa_track13"]
df["diff_activity_probability_sure"]=df["isa_track7"]-df["isa_track15"]
df["diff_cage_sure"]=df["isa_track1"]-df["isa_track7"]

# group by protein, and get the mean of the differences
df_grouped=df[['diff_activity_probability_cage', 
               'diff_activity_probability_dhs', 
               'diff_activity_probability_starr', 
               'diff_activity_probability_sure', 
               'diff_cage_sure', 
               'protein']].groupby('protein').mean()


df_med=pd.read_csv(f"2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)
tfs=df_med["gene"].unique().tolist()


df_grouped["in_med"]=df_grouped.index.isin(tfs).astype(int)


# kde plot
for col in df_grouped.columns[:-1]:
    plt.figure()
    sns.kdeplot(data=df_grouped[df_grouped["in_med"]==1][col],label="Mediator interactors")
    sns.kdeplot(data=df_grouped[df_grouped["in_med"]==0][col],label="Not interactors")
    # add p value of mann whitney u test
    u,p=mannwhitneyu(df_grouped[df_grouped["in_med"]==1][col],df_grouped[df_grouped["in_med"]==0][col])
    plt.text(0.5,0.5,f"p={p:.3f}",transform=plt.gca().transAxes)
    plt.title(col)
    plt.legend()
    plt.savefig(f"kde_{col}_{re}.png")
    plt.close()





#---------------------------------------
# is the difference related to histone modifications?
#---------------------------------------


df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{re}_k562.csv")

df["diff_activity_probability_cage"]=df["isa_track1"]-df["isa_track9"]
df["diff_activity_probability_dhs"]=df["isa_track3"]-df["isa_track11"]
df["diff_activity_probability_starr"]=df["isa_track5"]-df["isa_track13"]
df["diff_activity_probability_sure"]=df["isa_track7"]-df["isa_track15"]
df["diff_cage_sure"]=df["isa_track1"]-df["isa_track7"]


df_histone=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/{re}_k562.csv")
df_histone["end"]=df_histone["end"]-1
# add region in format chr:start-end
df_histone["region"]=df_histone["chrom"].astype(str)+":"+df_histone["start"].astype(str)+"-"+df_histone["end"].astype(str)
# merge with df by region
df=pd.merge(df,df_histone,on="region",how="inner")


# get columns starting with "log_signal"
histone_modifications=[col for col in df_histone.columns if "log_signal" in col]
histone_modifications=[col.replace("log_signal_","") for col in histone_modifications]


# calculate pearson correlation between differences in activity probability and histone modifications



def get_corr_all_histones(df,diff,histone_modifications):
    """
    Return a dataframe with one row. 
    Row name is df.protein, 
    columns are histone modifications, 
    values are pearson correlation coefficients if p is significant, 0 otherwise.
    """
    rs=[]
    for col in histone_modifications:
        r,p=pearsonr(df[diff],df["log_signal_"+col])
        rs.append(r if p<0.05 else 0)
    return pd.DataFrame([rs],columns=histone_modifications,index=[df.protein.unique()[0]])



def get_corr_all_tfs_all_histones(df,diff,histone_modifications):
    df_res=pd.DataFrame()
    for tf in df.protein.unique():
        df_sub=df[df["protein"]==tf].reset_index(drop=True)
        num_samples=df_sub.shape[0]
        df_temp=get_corr_all_histones(df_sub,diff,histone_modifications)
        df_temp["num_samples"]=num_samples
        df_res=pd.concat([df_res,get_corr_all_histones(df_sub,diff,histone_modifications)])
    df_res.to_csv(f"correlation_{diff}_{re}.csv")


def plot_heatmap(diff):
    df_res=pd.read_csv(f"correlation_{diff}_{re}.csv",index_col=0)
    # remove rows with all zeros
    df_res["row_abs_sum"]=df_res.abs().sum(axis=1)
    df_res=df_res[df_res["row_abs_sum"]>0].drop(columns="row_abs_sum")
    # clustermap
    sns.clustermap(df_res,cmap="coolwarm",center=0,figsize=(10,100))
    plt.savefig(f"heatmap_{diff}_{re}.pdf")
    plt.close()




for diff in ['diff_activity_probability_cage', 
               'diff_activity_probability_dhs', 
               'diff_activity_probability_starr', 
               'diff_activity_probability_sure', 
               'diff_cage_sure']:
    logger.info(diff)
    get_corr_all_tfs_all_histones(df,diff,histone_modifications)
    plot_heatmap(diff)

