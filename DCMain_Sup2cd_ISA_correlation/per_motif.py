import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
#----------------
# Helper functions
#----------------

def read_file(file_name):
    motif_df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{file_name}.csv",index_col=0)
    for track_num in range(8):
        motif_df[f"max_gradxinp_{track_num}"]=motif_df[f'gradxinp_{track_num}'].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        motif_df[f"max_ism_{track_num}"]=motif_df[f'ism_{track_num}'].apply(lambda x: max([float(i) for i in str(x).split(":")]))
        motif_df.drop(f'gradxinp_{track_num}',axis=1,inplace=True)
    # select only columns starting with "isa" or "max_gradxinp"
    cols_isa=[f"isa_track{track_num}" for track_num in range(8)]
    cols_max_gradxinp=[col for col in motif_df.columns if col.startswith("max_gradxinp")]
    cols_max_ism=[col for col in motif_df.columns if col.startswith("max_ism")]
    isa_flat=motif_df[cols_isa].values.flatten()
    max_gradxinp_flat=motif_df[cols_max_gradxinp].values.flatten()
    max_ism_flat=motif_df[cols_max_ism].values.flatten()
    df_res=pd.DataFrame({"isa":isa_flat,"max_gradxinp":max_gradxinp_flat,"max_ism":max_ism_flat})
    # shuffle rows
    return df_res.sample(frac=1).reset_index(drop=True)

def get_corr(df):
    corr_list_isa_max_gradxinp=[]
    corr_list_ism_max_gradxinp=[]
    corr_list_isa_ism=[]
    step=10000
    for i in range(0,df.shape[0],step):
        if i+step>df.shape[0]:
            step=df.shape[0]-i
        # calculate pearson correlation
        corr_isa_max_gradxinp,_=pearsonr(df["isa"][i:i+step],df["max_gradxinp"][i:i+step])
        corr_ism_max_gradxinp,_=pearsonr(df["max_ism"][i:i+step],df["max_gradxinp"][i:i+step])
        corr_isa_ism,_=pearsonr(df["isa"][i:i+step],df["max_ism"][i:i+step])
        corr_list_isa_max_gradxinp.append(corr_isa_max_gradxinp)
        corr_list_ism_max_gradxinp.append(corr_ism_max_gradxinp)
        corr_list_isa_ism.append(corr_isa_ism)
    return pd.DataFrame({"corr_isa_max_gradxinp":corr_list_isa_max_gradxinp,
                         "corr_ism_max_gradxinp":corr_list_ism_max_gradxinp,
                         "corr_isa_ism":corr_list_isa_ism})

#----------------
# Data generation
#----------------
df_promoters_hepg2=read_file("promoters_hepg2") # 202713
df_enhancers_hepg2=read_file("enhancers_hepg2") # 230086
df_promoters_k562=read_file("promoters_k562") # 271749
df_enhancers_k562=read_file("enhancers_k562") # 211824

df=pd.concat([df_promoters_hepg2,df_enhancers_hepg2,df_promoters_k562,df_enhancers_k562],axis=0)
# shuffle row
df=df.sample(frac=1).reset_index(drop=True)
df_corr=get_corr(df)


# save to csv
df_corr.to_csv("corr_per_motif.csv",index=False)

#---------------
# Plot
#---------------

def plot(title,col):
    df_corr=pd.read_csv("corr_per_motif.csv")
    # plot the correlation
    plt.figure(figsize=(3,2.5))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    sns.kdeplot(df_corr[col],linewidth=1)
    plt.title(title,fontsize=7)
    plt.xlabel("Pearson correlation",fontsize=7)
    plt.ylabel("Density",fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.xlim(0,1)
    plt.tight_layout()
    plt.savefig(f"{col}.pdf")
    plt.close()



plot("Motif level ISA vs max(gradxinp)","corr_isa_max_gradxinp")
plot("Motif level ISM vs max(gradxinp)","corr_ism_max_gradxinp")
plot("Motif level ISA vs ISM","corr_isa_ism")


# nohup python3 per_motif.py > per_motif.log &