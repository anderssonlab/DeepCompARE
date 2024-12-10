import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

file_path="/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/"

#-------------------------------------------------------------------------------
# Method 1: concatenate everything into 1d array, then calculate overall correlation
#-------------------------------------------------------------------------------

def read_files(file_suffix):
    df_isa=pd.read_csv(f"{file_path}isa_{file_suffix}.csv",index_col=0,header=None)
    df_isa["seq_idx"]=df_isa.index.str.split("_").str[0]
    df_ism=pd.read_csv(f"{file_path}ism_{file_suffix}.csv",index_col=0,header=None)
    df_ism["seq_idx"]=df_ism.index.str.split("_").str[0]
    df_gradxinp=pd.read_csv(f"{file_path}gradxinp_{file_suffix}.csv",index_col=0,header=None)
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


# randomly sample 10 sequences
def get_corr(df_ism,df_isa,df_gradxinp,max_seq_index, num_samples):
    sample_idx=np.random.choice(max_seq_index,num_samples,replace=False)
    sample_idx=[f"Seq{idx}" for idx in sample_idx]
    df_ism_sample=df_ism[df_ism["seq_idx"].isin(sample_idx)]
    df_isa_sample=df_isa[df_isa["seq_idx"].isin(sample_idx)]
    df_gradxinp_sample=df_gradxinp[df_gradxinp["seq_idx"].isin(sample_idx)]
    #
    assert all(df_ism_sample["seq_idx"]==df_isa_sample["seq_idx"])
    assert all(df_ism_sample["seq_idx"]==df_gradxinp_sample["seq_idx"])
    #
    ism_sample=df_ism_sample.iloc[:,:-1].values.flatten()
    isa_sample=df_isa_sample.iloc[:,:-1].values.flatten()
    gradxinp_sample=df_gradxinp_sample.iloc[:,:-1].values.flatten()
    #
    corr1,_=pearsonr(ism_sample,isa_sample)
    corr2,_=pearsonr(ism_sample,gradxinp_sample)
    corr3,_=pearsonr(isa_sample,gradxinp_sample)
    return corr1,corr2,corr3



def get_corr_all(df_ism,df_isa,df_gradxinp,max_seq_index,num_samples):
    corr1_list=[]
    corr2_list=[]
    corr3_list=[]
    for i in range(100):
        corr1,corr2,corr3=get_corr(df_ism,df_isa,df_gradxinp,max_seq_index,num_samples)
        corr1_list.append(corr1)
        corr2_list.append(corr2)
        corr3_list.append(corr3)
    df_res=pd.DataFrame({"corr_ism_isa":corr1_list,"corr_ism_gradxinp":corr2_list,"corr_isa_gradxinp":corr3_list})
    return df_res


def analyze(file_suffix,num_samples):
    df_ism,df_isa,df_gradxinp,max_seq_index=read_files(file_suffix)
    df_res=get_corr_all(df_ism,df_isa,df_gradxinp,max_seq_index,num_samples)
    df_res.to_csv(f"corr_base_{file_suffix}.csv")


analyze("promoters_hepg2",1000)
analyze("enhancers_hepg2",1000)
analyze("promoters_k562",1000)
analyze("enhancers_k562",1000)




# Plot


df_promoters_hepg2=pd.read_csv("corr_base_promoters_hepg2.csv",index_col=0)
df_promoters_hepg2["dataset"]="promoters_hepg2"
df_enhancers_hepg2=pd.read_csv("corr_base_enhancers_hepg2.csv",index_col=0)
df_enhancers_hepg2["dataset"]="enhancers_hepg2"
df_promoters_k562=pd.read_csv("corr_base_promoters_k562.csv",index_col=0)
df_promoters_k562["dataset"]="promoters_k562"
df_enhancers_k562=pd.read_csv("corr_base_enhancers_k562.csv",index_col=0)
df_enhancers_k562["dataset"]="enhancers_k562"


df=pd.concat([df_promoters_hepg2,df_enhancers_hepg2,df_promoters_k562,df_enhancers_k562])

sns.violinplot(data=df,x="dataset",y="corr_ism_isa")
plt.xlabel("Dataset")
plt.xticks(rotation=45)
plt.ylabel("Pearson correlation between ISA and ISM")
plt.tight_layout()
plt.savefig("violin_corr_ism_isa.pdf")
plt.close()



sns.violinplot(data=df,x="dataset",y="corr_isa_gradxinp")
plt.xlabel("Dataset")
plt.xticks(rotation=45)
plt.ylabel("Pearson correlation between ISA and gradxinp")
plt.tight_layout()
plt.savefig("violin_corr_isa_gradxinp.pdf")
plt.close()









#-------------------------------------------------------------------------------
# Method 2: calculate correlation in one sequence
#-------------------------------------------------------------------------------


file_suffix="promoters_hepg2"
df_ism,df_isa,df_gradxinp,max_seq_index=read_files(file_suffix)

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
    row_gradxinp=df_gradxinp.iloc[row_idx,:-1].values
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
plt.xlabel("Per sequence Pearson correlation between ISM and gradxinp")
plt.ylabel("Frequency")
plt.savefig("hist_per_sequence_corr_ism_gradxinp.pdf")
plt.close()

sns.histplot(corr_list_isa_gradxinp)
plt.xlabel("Per sequence Pearson correlation between ISA and gradxinp")
plt.ylabel("Frequency")
plt.savefig("hist_per_sequence_corr_isa_gradxinp.pdf")
plt.close()

















# nohup python3 per_base_importance_sequence_structure.py > per_base_importance_sequence_structure.log &