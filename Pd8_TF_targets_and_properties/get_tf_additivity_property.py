"""
From /isdata/alab/people/pcr980/DeepCompare/Motif_cooperitivity/df_sup_sub.csv
get name list of sub and super TFs
"""

import pandas as pd
import numpy as np

def split_dimer(tf_list):
    res_list=[]
    for tf in tf_list:
        if "::" not in tf:
            res_list.append(tf)
        else:
            res_list.extend(tf.split("::"))
    return res_list

def write():
    with open("sub_tfs.txt","w") as f:
        f.write("\n".join(sub_tfs))
    with open("super_tfs.txt","w") as f:
        f.write("\n".join(super_tfs))

# read in additivity profile
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Motif_cooperitivity/df_sup_sub.csv")
df=df[(df["sig_super_count"]+df["sig_sub_count"])>10].reset_index(drop=True)
df["track_num"]=df["track_num"].map({0: 'cage', 1: 'cage', 2: 'dhs', 3: 'dhs', 4: 'starr', 5: 'starr', 6: 'sure', 7: 'sure'})

# write TF property for promoters/enhancers hepg2/k562
df["super_add"]=df["super_sub_ratio"] > 0.5
df["sub_add"]=df["super_sub_ratio"] < -0.5
df=df.groupby(["dataset","protein"]).agg({"super_add":"sum","sub_add":"sum"}).reset_index()
df.rename(columns={"super_add":"super_profile_count","sub_add":"sub_profile_count"},inplace=True)
df["tf_property"]=np.where(df["super_profile_count"]==4,"super",np.where(df["sub_profile_count"]>=3,"sub","unknown"))
df_additivity=df.pivot(index="protein",columns="dataset",values="tf_property").reset_index()
df_additivity.to_csv("tf_additivity_property.csv",index=False)


# classify tf as sub, super
# define column "TF_type": 
# sub: count_sub>1, count_super==0
# super: count_super>0, count_sub==0

df_additivity["count_sub"]=df_additivity.apply(lambda x: x.isin(["sub"]).sum(),axis=1)
df_additivity["count_super"]=df_additivity.apply(lambda x: x.isin(["super"]).sum(),axis=1)
df_additivity["TF_type"]="unknown"
df_additivity.loc[(df_additivity["count_sub"]>1) & (df_additivity["count_super"]==0),"TF_type"]="sub"
df_additivity.loc[(df_additivity["count_super"]>0) & (df_additivity["count_sub"]==0),"TF_type"]="super"

sub_tfs=df_additivity[df_additivity["TF_type"]=="sub"]["protein"].to_list()
super_tfs=df_additivity[df_additivity["TF_type"]=="super"]["protein"].to_list()
super_tfs=set(split_dimer(super_tfs))
write()


