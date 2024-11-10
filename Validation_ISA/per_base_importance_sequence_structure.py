import pandas as pd
import numpy as np
from scipy.stats import pearsonr

file_suffix="enhancers_hepg2"
file_path="/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/"

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


# randomly sample 10 sequences

sample_idx=np.random.choice(max_seq_index1,1000,replace=False)
sample_idx=[f"Seq{idx}" for idx in sample_idx]
df_ism_sample=df_ism[df_ism["seq_idx"].isin(sample_idx)]
df_isa_sample=df_isa[df_isa["seq_idx"].isin(sample_idx)]
df_gradxinp_sample=df_gradxinp[df_gradxinp["seq_idx"].isin(sample_idx)]

assert all(df_ism_sample["seq_idx"]==df_isa_sample["seq_idx"])
assert all(df_ism_sample["seq_idx"]==df_gradxinp_sample["seq_idx"])

ism_sample=df_ism_sample.iloc[:,:-1].values.flatten()
isa_sample=df_isa_sample.iloc[:,:-1].values.flatten()
gradxinp_sample=df_gradxinp_sample.iloc[:,:-1].values.flatten()
pearsonr(ism_sample,isa_sample)
pearsonr(ism_sample,gradxinp_sample)
pearsonr(isa_sample,gradxinp_sample)