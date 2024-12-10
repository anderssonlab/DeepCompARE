import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from loguru import logger

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")

#----------------
# Helper functions and tools
#----------------

def get_corr(motif_df):
    corr_list=[]
    for track_num in range(8):
        corr,_=pearsonr(motif_df[f'isa_track{track_num}'],motif_df[f'max_gradxinp_{track_num}'])
        corr_list.append(corr)
    return corr_list


file_name="promoters_hepg2"
def get_corr_repeats(file_name):
    motif_df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_name}.csv",index_col=0)
    for track_num in range(8):
        motif_df[f"max_gradxinp_{track_num}"]=motif_df[f'gradxinp_{track_num}'].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        motif_df.drop(f'gradxinp_{track_num}',axis=1,inplace=True)
    corr_list=[]
    for i in range(100):
        df_sample=motif_df.iloc[:1000,].reset_index(drop=True)
        # df_sample=motif_df.sample(1000).reset_index(drop=True)
        corr=get_corr(df_sample)
        corr_list.append(corr)
    return corr_list




# plt.figure(figsize=(8,6))
# track_num=0
# sns.scatterplot(data=df_sample,x=f'isa_track{track_num}',y=f'max_gradxinp_{track_num}')
# plt.xlabel("Motif ISA")
# plt.ylabel("Max(gradxinp)")
# plt.title(f"Track {track_num}")
# plt.savefig(f"scatter_continuous.pdf")
# plt.close()


#---------
# Data generation
#----------------

corr_list=[]
for file_name in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
    corr_list+=get_corr_repeats(file_name)

corr_df=pd.DataFrame(corr_list,columns=[f"track_{i}" for i in range(16)])
corr_df.to_csv("corr_motif_isa_vs_max_gradxinp.csv",index=False)


#---------------
# Plot
#---------------


# df=pd.read_csv("corr_motif_isa_vs_max_gradxinp.csv")
# df=df.melt(var_name="track",value_name="corr")
# # kde plot, hue=track
# plt.figure(figsize=(8,6))
# sns.kdeplot(data=df,x="corr",hue="track")
# plt.xlabel("Pearson correlation between motif ISA and max(gradxinp)")
# plt.savefig("corr_motif_isa_vs_max_gradxinp.pdf")
# plt.close()






# nohup python3 motif_isa_vs_max_gradxinp.py > motif_isa_vs_max_gradxinp.log &