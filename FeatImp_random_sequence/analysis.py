
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
    df_isa=pd.read_csv(f"Pd1_random_seqs_and_feat_imp/df_isa.csv",index_col=0,header=None)
    df_isa["seq_idx"]=df_isa.index.str.split("_").str[0]
    df_ism=pd.read_csv(f"Pd1_random_seqs_and_feat_imp/df_ism.csv",index_col=0,header=None)
    df_ism["seq_idx"]=df_ism.index.str.split("_").str[0]
    df_gradxinp=pd.read_csv(f"Pd1_random_seqs_and_feat_imp/df_gradxinp.csv",index_col=0,header=None)
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
    df_ism.reset_index(drop=True,inplace=True)
    df_isa.reset_index(drop=True,inplace=True)  
    df_gradxinp.reset_index(drop=True,inplace=True)
    return df_ism,df_isa,df_gradxinp,max_seq_index1




df_ism,df_isa,df_gradxinp,max_seq_index=read_files()

#---------------------------
# Per basexsequence correlation
#---------------------------
# overall correlation
ism_flat=df_ism.iloc[:,:-1].values.flatten()
isa_flat=df_isa.iloc[:,:-1].values.flatten()
gradxinp_flat=df_gradxinp.iloc[:,:-1].values.flatten()
corr_ism_isa, _=pearsonr(ism_flat,isa_flat) # 0.913
corr_ism_gradxinp, _=pearsonr(ism_flat,gradxinp_flat) # -0.0001
corr_isa_gradxinp, _=pearsonr(isa_flat,gradxinp_flat) # 0.003
corr_ism_abs_gradxinp, _=pearsonr(ism_flat,np.abs(gradxinp_flat)) # -0.003
corr_isa_abs_gradxinp, _=pearsonr(isa_flat,np.abs(gradxinp_flat)) # 0.01






#---------------------------
# Per sequence correlation
#---------------------------

corr_list_ism_isa=[]
corr_list_ism_gradxinp=[]
corr_list_isa_gradxinp=[]
corr_list_ism_abs_gradxinp=[]
corr_list_isa_abs_gradxinp=[]

for row_idx in range(df_gradxinp.shape[0]):
    row_ism=df_ism.iloc[row_idx,:-1].values
    row_isa=df_isa.iloc[row_idx,:-1].values
    row_gradxinp=df_gradxinp.iloc[row_idx,:-1].values
    row_gradxinp_abs=np.abs(row_gradxinp)
    corr1,_=pearsonr(row_ism,row_isa)
    corr2,_=pearsonr(row_ism,row_gradxinp)
    corr3,_=pearsonr(row_isa,row_gradxinp)
    corr4,_=pearsonr(row_ism,row_gradxinp_abs)
    corr5,_=pearsonr(row_isa,row_gradxinp_abs)
    corr_list_ism_isa.append(corr1)
    corr_list_ism_gradxinp.append(corr2)
    corr_list_isa_gradxinp.append(corr3)
    corr_list_ism_abs_gradxinp.append(corr4)
    corr_list_isa_abs_gradxinp.append(corr5)


def plot_kde(corr_list, title, filename):
    plt.figure(figsize=(3, 2.5))
    sns.kdeplot(corr_list,linewidth=1)
    plt.xlabel("Pearson correlation",fontsize=7)
    plt.ylabel("Frequency",fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.title(title,fontsize=7)
    plt.xlim(0,1)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


plot_kde(corr_list_ism_isa, "ISM vs ISA", "kde_ism_isa.pdf")
plot_kde(corr_list_ism_gradxinp, "ISM vs gradxinp", "kde_ism_gradxinp.pdf")
plot_kde(corr_list_isa_gradxinp, "ISA vs gradxinp", "kde_isa_gradxinp.pdf")
plot_kde(corr_list_ism_abs_gradxinp, "ISM vs abs(gradxinp)", "kde_ism_abs_gradxinp.pdf")
plot_kde(corr_list_isa_abs_gradxinp, "ISA vs abs(gradxinp)", "kde_isa_abs_gradxinp.pdf")

