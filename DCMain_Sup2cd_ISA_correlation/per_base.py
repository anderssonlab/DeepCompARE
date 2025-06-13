import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




file_path="/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/"


# generate shuffling index to shuffle rows 
def shuffle(df_ism,df_isa,df_gradxinp):
    shuffling_index=np.random.choice(df_ism.shape[0],df_ism.shape[0],replace=False)
    df_ism=df_ism.iloc[shuffling_index,:]
    df_isa=df_isa.iloc[shuffling_index,:]
    df_gradxinp=df_gradxinp.iloc[shuffling_index,:]
    return df_ism,df_isa,df_gradxinp


def get_corr_per_base():
    corr_list_ism_isa=[]
    corr_list_ism_gradxinp=[]
    corr_list_isa_gradxinp=[]
    step=10000
    for i in range(0,df_ism.shape[0],step):
        logger.info(f"Processing {i} to {i+step} out of {df_ism.shape[0]}")
        if i+step>df_ism.shape[0]:
            break
        row_ism=df_ism.iloc[i:i+step,:-1].values.flatten()
        row_isa=df_isa.iloc[i:i+step,:-1].values.flatten()
        row_gradxinp=df_gradxinp.iloc[i:i+step,:-1].values.flatten()
        corr1,_=pearsonr(row_ism,row_isa)
        corr2,_=pearsonr(row_ism,row_gradxinp)
        corr3,_=pearsonr(row_isa,row_gradxinp)
        corr_list_ism_isa+=[corr1]
        corr_list_ism_gradxinp+=[corr2]
        corr_list_isa_gradxinp+=[corr3]
    return corr_list_ism_isa,corr_list_ism_gradxinp,corr_list_isa_gradxinp


corr_list_ism_isa=[]
corr_list_ism_gradxinp=[]
corr_list_isa_gradxinp=[]

for file_name in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]: # 16775+15238+12614+10450 
    logger.info(f"Processing {file_name}")
    df_ism=pd.read_csv(f"{file_path}/ism_{file_name}.csv",index_col=0).reset_index(drop=True).round(3)
    df_isa=pd.read_csv(f"{file_path}/isa_{file_name}.csv",index_col=0).reset_index(drop=True).round(3)
    df_gradxinp=pd.read_csv(f"{file_path}/gradxinp_{file_name}.csv",index_col=0).reset_index(drop=True).round(3)
    # shuffle
    df_ism,df_isa,df_gradxinp=shuffle(df_ism,df_isa,df_gradxinp)
    # get correlation
    corr_list_ism_isa_temp,corr_list_ism_gradxinp_temp,corr_list_isa_gradxinp_temp=get_corr_per_base()
    corr_list_ism_isa+=corr_list_ism_isa_temp
    corr_list_ism_gradxinp+=corr_list_ism_gradxinp_temp
    corr_list_isa_gradxinp+=corr_list_isa_gradxinp_temp

# output the correlation
df_corr=pd.DataFrame({"corr_ism_isa":corr_list_ism_isa,"corr_ism_gradxinp":corr_list_ism_gradxinp,"corr_isa_gradxinp":corr_list_isa_gradxinp})
df_corr.to_csv("corr_per_base.csv",index=False)


df_corr=pd.read_csv("corr_per_base.csv")

def plot(title,col):
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

plot("Base level ISA v.s. ISM","corr_ism_isa")
plot("Base level ISA v.s. gradxinp","corr_isa_gradxinp")
plot("Base level ISM v.s. gradxinp","corr_ism_gradxinp")







# nohup python3 per_base.py > per_base.log &