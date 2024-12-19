import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


cell_line="k562"

df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
df=df[df["c_sum"]>1].reset_index(drop=True)

# get tf type
tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt",header=None).iloc[:,0].tolist()
tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt",header=None).iloc[:,0].tolist()


df["protein2_cooperativity"]="other"
df.loc[df["protein2"].isin(tfs_codependent),"protein2_cooperativity"]="codependent"
df.loc[df["protein2"].isin(tfs_redundant),"protein2_cooperativity"]="redundant"





df_sub=df[df["protein2"]=="BACH1"].reset_index(drop=True)
df_sub=df_sub[["protein2","c_redundancy","c_codependency"]]


# long format
df_sub=df_sub.melt(id_vars=["protein2"],value_vars=["c_redundancy","c_codependency"],var_name="cooperativity_type",value_name="cooperativity_index")
# remove columns cooperativity_type
df_sub=df_sub.drop(columns=["cooperativity_type"])
# get the rank of the cooperativity index
df_sub["rank"]=df_sub["cooperativity_index"].rank(ascending=False)


# scatterplot rank v.s. cooperativity index

sns.scatterplot(data=df_sub,x="rank",y="cooperativity_index")
plt.xlabel("Rank",fontsize=15)
plt.ylabel("Cooperativity Index",fontsize=15)
plt.savefig(f"rank_cooperativity_index.png")
plt.close()
