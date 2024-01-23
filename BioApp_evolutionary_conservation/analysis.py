import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import extract_numbers




#-------------------------------------------------------------------------------------
# Analysis 1: Pearson correlation between feature importance and conservation
#-------------------------------------------------------------------------------------
track_info={0: "CAGE_HepG2",
            1: "CAGE_K562",
            2: "DHS_HepG2",
            3: "DHS_K562",
            4: "STARR_HepG2",
            5: "STARR_K562",
            6: "SuRE_HepG2",
            7: "SuRE_K562"}

def read_featimp(featimp_file,track_num):
    """
    Read featimp from featimp_file, subset by track_num
    Featimp is either gradxinp or ism
    """
    featimp_df=pd.read_csv(featimp_file,index_col=0)
    # Given that indices are composed of "SeqX_TrackY", we can subset to contain only "_Track{track_num}"
    featimp_df=featimp_df[featimp_df.index.str.contains(f"_Track{track_num}$")]
    # reorder rows by number 
    sorted_indices = sorted(featimp_df.index, key=lambda x: extract_numbers(x))
    featimp_df=featimp_df.reindex(sorted_indices)
    assert np.all(featimp_df.index==f'Seq{i}_Track{track_num}' for i in range(len(featimp_df)))
    return featimp_df


# def compute_correlation(feat_imp_file,cons_file,track_num,absolute_importance):
#     """
#     Args:
#         feat_imp_file: file path to feature importance file
#         cons_file: file path to conservation file
#         track_num: track number to use
#     Return:
#         corr: correlation between feature importance and conservation
#         pval: p-value
#     """
#     df_feat_imp=read_featimp(feat_imp_file,track_num=track_num)
#     df_cons=pd.read_csv(cons_file,header=None)
#     feat_imp=df_feat_imp.values.reshape(-1,1).squeeze()
#     if absolute_importance:
#         feat_imp=np.abs(feat_imp)
#     cons=df_cons.values.reshape(-1,1).squeeze()
#     assert len(feat_imp)==len(cons)
#     mask_feat_imp_nan = np.isnan(feat_imp)
#     mask_feat_imp_inf = np.isinf(feat_imp)
#     mask_cons_nan = np.isnan(cons)
#     mask_cons_inf = np.isinf(cons)
#     mask_either = mask_feat_imp_nan | mask_feat_imp_inf | mask_cons_nan | mask_cons_inf
#     corr,pval=pearsonr(feat_imp[~mask_either],cons[~mask_either])
#     return corr,pval



# corrs=[]
# pvals=[]
# tracks=[]
# absolute_importances=[]
# conservation_files=[]
# for i in range(8):
#     for absolute_importance in [True,False]:
#         for conservation_file in ["phyloP","phastCons"]:
#             logger.info(f"Track: {track_info[i]}, absolute_importance: {absolute_importance}, conservation_file: {conservation_file}")
#             corr,pval=compute_correlation(
#                 f"/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/ism_{track_info[i]}.csv", # or graxinp_{track_info[i]}.csv
#                 f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/{conservation_file}_{track_info[i]}.csv",
#                 track_num=i,
#                 absolute_importance=absolute_importance)
#             corrs.append(corr)
#             pvals.append(pval)
#             tracks.append(track_info[i])
#             absolute_importances.append(absolute_importance)
#             conservation_files.append(conservation_file)
# pd.DataFrame({"track":tracks,         
#               "absolute_importance":absolute_importances,
#               "conservation_file":conservation_files,
#               "correlation":corrs,
#               "p-value":pvals}).to_csv("correlation_ismn_vs_conservation.csv",index=False)
    



#-------------------------------------------------------  
# Analysis 2: plot correlations
#-------------------------------------------------------  

# df_corr=pd.read_csv("correlation_ismn_vs_conservation.csv")
# df_sub=df_corr[df_corr["absolute_importance"]==True]
# fig=plt.figure(figsize=(8,8))
# sns.barplot(x="track",y="correlation",hue="conservation_file",data=df_sub)
# plt.xticks(rotation=45)
# plt.title("ISM-N:PhastCons correlates more than phylop")
# plt.savefig("ismn_phastcon_vs_phylop.png")
# plt.close()

# df_sub=df_corr[df_corr["conservation_file"]=="phastCons"]
# fig=plt.figure(figsize=(8,8))
# sns.barplot(x="track",y="correlation",hue="absolute_importance",data=df_sub)
# plt.xticks(rotation=45)
# plt.title("phastCon: absolute ism correlates more than ism")
# plt.savefig("ismn_abs_phastcon.png")
# plt.close()



#-------------------------------------------------------  
# Analysis 3: bin ism values and compute mean conservation for each bin
#-------------------------------------------------------  

def calculate_mean_cons(conservation_file, track_num):
    df_feat_imp=read_featimp(f"/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/ism_{track_info[track_num]}.csv",track_num)
    df_cons=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/{conservation_file}_{track_info[track_num]}.csv",header=None)
    feat_imp=df_feat_imp.values.reshape(-1,1).squeeze()
    cons=df_cons.values.reshape(-1,1).squeeze()
    assert feat_imp.shape==cons.shape
    mask_feat_imp_nan = np.isnan(feat_imp)
    mask_feat_imp_inf = np.isinf(feat_imp)
    mask_cons_nan = np.isnan(cons)
    mask_cons_inf = np.isinf(cons)
    mask_either = mask_feat_imp_nan | mask_feat_imp_inf | mask_cons_nan | mask_cons_inf
    feat_imp=feat_imp[~mask_either]
    cons=cons[~mask_either]

    bins=[-np.inf,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,np.inf]
    feat_imp_bined = np.digitize(feat_imp, bins)
    data = pd.DataFrame({'Bin': feat_imp_bined, 'Cons': cons, 'ism':feat_imp})
    bin_labels = [f"{bins[i]:.2f} - {bins[i+1]:.2f}" for i in range(len(bins)-1)]
    data['Bin'] = data['Bin'].apply(lambda x: bin_labels[x-1])

    # calculate mean Cons for each bin
    data=data.groupby('Bin').median().reset_index()
    # Convert Bin to categorical
    data['Bin'] = pd.Categorical(data['Bin'], categories=bin_labels)
    data['conservation_file']=conservation_file
    data['track']=track_info[track_num]
    return data


data=pd.DataFrame()
for i in range(8):
    for conservation_file in ["phyloP","phastCons"]:
        logger.info(f"Track: {track_info[i]}, conservation_file: {conservation_file}")
        temp=calculate_mean_cons(conservation_file=conservation_file, track_num=i)
        data=pd.concat([data,temp])
data.to_csv("ismn_median_conservation.csv",index=False)



# plot mean Cons 
df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/ismn_median_conservation.csv")

# relevel df["Bin"] to be in order
bins=[-np.inf,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,np.inf]
bin_labels = [f"{bins[i]:.2f} - {bins[i+1]:.2f}" for i in range(len(bins)-1)]
df['Bin'] = pd.Categorical(df['Bin'], categories=bin_labels)

# scatter plot, x=bin, y=mean cons, color=track shape=conservation_file
fig=plt.figure(figsize=(10,8))
sns.scatterplot(x="Bin",y="Cons",hue="track",style="conservation_file",data=df)
plt.ylabel("Median conservation per bin")
plt.xticks(rotation=45) 
# put legend outside of plot, right next to plot
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0)
plt.title("Median conservation for each bin of ISM")
plt.savefig("ismn_bined_median_conservation.png")
plt.close()