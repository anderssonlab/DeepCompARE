import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import read_featimp



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
# feat_ipm_files=[]
# for i in range(8):
#     for absolute_importance in [True,False]:
#         for conservation_file in ["phyloP447way","phyloP447wayLRT"]:
#             for feat_ipm_file in ["ism","gradxinp"]:
#                 corr,pval=compute_correlation(
#                     f"/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/{feat_ipm_file}_{track_info[i]}.csv", 
#                     f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/{conservation_file}_{track_info[i]}.csv",
#                     track_num=i,
#                     absolute_importance=absolute_importance)
#                 logger.info(f"Track: {track_info[i]},conservation_file: {conservation_file}, feat_imp_file: {feat_ipm_file}, absolute_importance: {absolute_importance},Correlation: {corr},p-value: {pval}")
#                 corrs.append(corr)
#                 pvals.append(pval)
#                 tracks.append(track_info[i])
#                 absolute_importances.append(absolute_importance)
#                 conservation_files.append(conservation_file)
#                 feat_ipm_files.append(feat_ipm_file)
                
# pd.DataFrame({"track":tracks,         
#               "absolute_importance":absolute_importances,
#               "conservation_file":conservation_files,
#               "feat_ipm_file":feat_ipm_files,
#               "correlation":corrs,
#               "p-value":pvals}).to_csv("correlation_featimp_vs_conservation.csv",index=False)

# logger.info("Done with analysis 1.")



#-------------------------------------------------------  
# Analysis 2: plot correlations
#-------------------------------------------------------  

# df_corr=pd.read_csv("correlation_featimp_vs_conservation.csv")
# df_sub=df_corr[df_corr["absolute_importance"]==False]

# fig=plt.figure(figsize=(8,8))
# sns.scatterplot(x="track",y="correlation",hue="conservation_file",style="feat_ipm_file",data=df_sub)
# plt.xticks(rotation=45)
# plt.title("Feature importance vs conservation")
# plt.savefig("correlation_featimp_vs_conservation_no_abs.pdf")
# plt.close()



#-------------------------------------------------------  
# Analysis 3: bin ism values and compute mean conservation for each bin
#-------------------------------------------------------  

# use ism by default
def get_bin_df(conservation_file, track_num):
    df_ism=read_featimp(f"/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/ism_{track_info[track_num]}.csv",track_num)
    df_cons=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/{conservation_file}_{track_info[track_num]}.csv",header=None)
    ism=df_ism.values.reshape(-1,1).squeeze()
    cons=df_cons.values.reshape(-1,1).squeeze()
    assert ism.shape==cons.shape
    mask_ism_nan = np.isnan(ism)
    mask_ism_inf = np.isinf(ism)
    mask_cons_nan = np.isnan(cons)
    mask_cons_inf = np.isinf(cons)
    mask_either = mask_ism_nan | mask_ism_inf | mask_cons_nan | mask_cons_inf
    ism=ism[~mask_either]
    cons=cons[~mask_either]
    
    # set bins
    bins=[-np.inf,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,np.inf]
    ism_bined = np.digitize(ism, bins)
    data = pd.DataFrame({'Bin': ism_bined, 'Cons': cons, 'ism':ism})
    bin_labels = [f"{bins[i]:.2f} - {bins[i+1]:.2f}" for i in range(len(bins)-1)]
    data['Bin'] = data['Bin'].apply(lambda x: bin_labels[x-1])
    data['Bin'] = pd.Categorical(data['Bin'], categories=bin_labels)
    return data


def calculate_mean_cons(conservation_file, track_num):
    data=get_bin_df(conservation_file, track_num)

    # calculate mean Cons for each bin
    data=data.groupby('Bin').mean().reset_index()
    # Convert Bin to categorical
    data['conservation_file']=conservation_file
    data['track']=track_info[track_num]
    return data

# calculate mean conservation for each bin
# data=pd.DataFrame()
# for i in range(8):
#     for conservation_file in ["phyloP447way","phyloP447wayLRT"]:
#         logger.info(f"Track: {track_info[i]}, conservation_file: {conservation_file}")
#         temp=calculate_mean_cons(conservation_file=conservation_file, track_num=i)
#         data=pd.concat([data,temp])
# data.to_csv("ism_bined_mean_conservation.csv",index=False)
# logger.info("Done with analysis 3.")



# plot mean Cons 
# df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/BioApp_evolutionary_conservation/ism_bined_mean_conservation.csv")

# # relevel df["Bin"] to be in order
# bins=[-np.inf,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,np.inf]
# bin_labels = [f"{bins[i]:.2f} - {bins[i+1]:.2f}" for i in range(len(bins)-1)]
# df['Bin'] = pd.Categorical(df['Bin'], categories=bin_labels)

# fig=plt.figure(figsize=(10,8))
# sns.scatterplot(x="Bin",y="Cons",hue="track",style="conservation_file",data=df)
# plt.ylabel("Mean conservation per bin")
# plt.xticks(rotation=45) 
# # put legend outside of plot, right next to plot
# plt.legend(bbox_to_anchor=(1, 1), loc=5, borderaxespad=2)
# plt.title("Mean conservation for each bin of ISM")
# plt.savefig("ismn_bined_mean_conservation.pdf")
# plt.close()




for conservation_file in ["phyloP447way","phyloP447wayLRT"]:
    for i in range(8):
        data=get_bin_df(conservation_file, i)
        sns.violinplot(x="Bin",y="Cons",data=data)
        plt.xticks(rotation=45) 
        plt.title(f"{conservation_file} for each bin of ISM {track_info[i]}")
        plt.savefig(f"violin_{conservation_file}_distribution_per_ism_bin_{track_info[i]}.pdf")
        plt.close()