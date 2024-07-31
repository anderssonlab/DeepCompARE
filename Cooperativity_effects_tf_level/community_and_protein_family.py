import pandas as pd
from communities.algorithms import louvain_method


# convert tf_pair to redundancy graph and codependent graph
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_pair_cooperativity_ratio_post_filter.csv")
# how many TF pairs have cooperativity_ratio between 0.3 and 0.7
df[(df["cooperativity_ratio"]>0.3) & (df["cooperativity_ratio"]<0.7)].shape[0] 

df=df.pivot(index="protein1",columns="protein2",values="cooperativity_ratio")
protein_names=df.index.tolist()
mat=(df>0.9).astype(int)
community_idx, _ = louvain_method(mat.values)
# convert protein name accordingt to protein_names
res=[]
for this_set in community_idx:
    this_set=list(this_set)
    res_set=[]
    for i in range(len(this_set)):
        res_set.append(protein_names[this_set[i]])
    res.append(res_set)

# retain only list with more than 1 element
res=[x for x in res if len(x)>1]

res[0]






#------------------------------------------------
# Do TF pair coming from same family have lower cooperativity ratio?
#------------------------------------------------

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_individual_effect_and_cooperativity/tf_pair_cooperativity_ratio_post_filter.csv")
df_family=pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_tf_family.csv")
df_family["ID"]=df_family["ID"].str.upper()

# are all df["protein2"] in df_family["ID"]?
# annotate protein1 and protein2 with family using merge
df=pd.merge(df,df_family,left_on="protein1",right_on="ID",how="left")
# drop columns 'AC', 'ID'
# change 'tf_family' to 'family_1'
df.drop(columns=['AC', 'ID'],inplace=True)
df.rename(columns={"tf_family":"family_1"},inplace=True)
df=pd.merge(df,df_family,left_on="protein2",right_on="ID",how="left")
df.drop(columns=['AC', 'ID'],inplace=True)
df.rename(columns={"tf_family":"family_2"},inplace=True)
df["same_family"]= (df["family_1"]==df["family_2"])

# sns kde plot
import seaborn as sns
import matplotlib.pyplot as plt
sns.kdeplot(data=df,x="cooperativity_ratio",hue="same_family",common_norm=False)
plt.savefig("Plots/cooperativity_ratio_same_family.pdf")
plt.close()

# calculate p value
from scipy.stats import mannwhitneyu
df_same_family=df[df["same_family"]==True]
df_diff_family=df[df["same_family"]==False]
mannwhitneyu(df_same_family["cooperativity_ratio"],df_diff_family["cooperativity_ratio"])
df_same_family["cooperativity_ratio"].median()
df_diff_family["cooperativity_ratio"].median()