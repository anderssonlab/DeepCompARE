
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python/")
from seq_ops import generate_random_seqs

seqs=generate_random_seqs(10000, 600)
pd.DataFrame(seqs).to_csv("random_seqs.csv", index=False, header=False) 



def read_files():
    df_isa=pd.read_csv(f"df_isa.csv",index_col=0,header=None)
    df_isa["seq_idx"]=df_isa.index.str.split("_").str[0]
    df_ism=pd.read_csv(f"df_ism.csv",index_col=0,header=None)
    df_ism["seq_idx"]=df_ism.index.str.split("_").str[0]
    df_gradxinp=pd.read_csv(f"df_gradxinp.csv",index_col=0,header=None)
    df_gradxinp["seq_idx"]=df_gradxinp.index.str.split("_").str[0]
    # get the max index of the sequence
    max_seq_index1=df_isa.iloc[-1,:].name.split("_")[0].strip("Seq")
    max_seq_index1=int(max_seq_index1)
    max_seq_index2=df_ism.iloc[-1,:].name.split("_")[0].strip("Seq")
    max_seq_index2=int(max_seq_index2)
    max_seq_index3=df_gradxinp.iloc[-1,:].name.split("_")[0].strip("Seq")
    max_seq_index3=int(max_seq_index3)
    assert max_seq_index1==max_seq_index2
    assert max_seq_index1==max_seq_index3
    return df_ism,df_isa,df_gradxinp,max_seq_index1




df_ism,df_isa,df_gradxinp,max_seq_index=read_files()

df_ism.reset_index(drop=True,inplace=True)
df_isa.reset_index(drop=True,inplace=True)
df_gradxinp.reset_index(drop=True,inplace=True)

sample_idx=np.random.choice(df_gradxinp.shape[0],1000,replace=False)


corr_list_ism_isa=[]
corr_list_ism_gradxinp=[]
corr_list_isa_gradxinp=[]


for row_idx in sample_idx:
    row_ism=df_ism.iloc[row_idx,:-1].values
    row_isa=df_isa.iloc[row_idx,:-1].values
    row_gradxinp=df_gradxinp.iloc[row_idx,:-1].abs().values
    corr1,_=pearsonr(row_ism,row_isa)
    corr2,_=pearsonr(row_ism,row_gradxinp)
    corr3,_=pearsonr(row_isa,row_gradxinp)
    corr_list_ism_isa.append(corr1)
    corr_list_ism_gradxinp.append(corr2)
    corr_list_isa_gradxinp.append(corr3)


sns.histplot(corr_list_ism_isa)
plt.xlabel("Per sequence Pearson correlation between ISM and ISA")
plt.ylabel("Frequency")
plt.savefig("hist_per_sequence_corr_ism_isa.pdf")
plt.close()

sns.histplot(corr_list_ism_gradxinp)
plt.xlabel("Per sequence Pearson correlation between ISM and |gradxinp|")
plt.ylabel("Frequency")
plt.savefig("hist_per_sequence_corr_ism_abs_gradxinp.pdf")
plt.close()

sns.histplot(corr_list_isa_gradxinp)
plt.xlabel("Per sequence Pearson correlation between ISA and |gradxinp|")
plt.ylabel("Frequency")
plt.savefig("hist_per_sequence_corr_isa_abs_gradxinp.pdf")
plt.close()



