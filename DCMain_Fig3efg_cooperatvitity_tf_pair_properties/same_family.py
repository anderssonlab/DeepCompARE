import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
# Do TF pair coming from same family have lower cooperativity ratio?

# TODO: change to merged
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_merged.csv")
df=df[df["c_sum"]>1].reset_index(drop=True)
df_family=pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_tf_family.csv")
df_family["ID"]=df_family["ID"].str.upper()

# are all df["protein2"] in df_family["ID"]?
# annotate protein1 and protein2 with family using merge
df=pd.merge(df,df_family,left_on="protein1",right_on="ID",how="inner")
df.drop(columns=['AC', 'ID'],inplace=True)
df.rename(columns={"tf_family":"family_1"},inplace=True)
df=pd.merge(df,df_family,left_on="protein2",right_on="ID",how="inner")
df.drop(columns=['AC', 'ID'],inplace=True)
df.rename(columns={"tf_family":"family_2"},inplace=True)
df["same_family"]= (df["family_1"]==df["family_2"])

# retain the "protein2" with at least 4 pair with same_family==True
df_filtered=pd.DataFrame()
for family,df_sub in df.groupby("family_2"):
    if df_sub["same_family"].sum()>=10:
        df_filtered=pd.concat([df_filtered,df_sub],axis=0)

# sns violin plot the distributon of cooperativity_index split by family
plt.figure(figsize=(10,10))
sns.violinplot(data=df_filtered, x='family_2', y='cooperativity_index',hue="same_family", split=True, cut=0, quantiles=[0.25, 0.5, 0.75],inner=None)
plt.xticks(rotation=90)
# make bottom margin larger
plt.gcf().subplots_adjust(bottom=0.7)
plt.savefig("same_family.pdf")
plt.close()





        
df_family["tf_family"].unique()
df_family[df_family["tf_family"]=="Three-zinc finger Kruppel-related"]
