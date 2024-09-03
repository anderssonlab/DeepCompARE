import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
# Do TF pair coming from same family have lower cooperativity ratio?

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_pair_cooperativity_ratio_post_filter_lenient.csv")
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
sns.kdeplot(data=df,x="cooperativity_ratio",hue="same_family",common_norm=False)
plt.savefig("Plots/cooperativity_ratio_same_family_lenient.pdf")
plt.close()

# calculate p value
df_same_family=df[df["same_family"]==True]
df_diff_family=df[df["same_family"]==False]
mannwhitneyu(df_same_family["cooperativity_ratio"],df_diff_family["cooperativity_ratio"])
df_same_family["cooperativity_ratio"].median()
df_diff_family["cooperativity_ratio"].median()