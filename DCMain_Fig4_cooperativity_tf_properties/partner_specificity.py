import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_merged.csv")
df=df[df["c_sum"]>1].reset_index(drop=True)
df["cooperativity"]=pd.cut(df["cooperativity_index"],bins=[-1,0.3,0.7,2],labels=["redundant","other","codependent"])

# group by protein2 and cooperativity, count the number of  protein1
df_coop=df.groupby(["protein2","cooperativity"]).agg({"protein1":pd.Series.nunique}).reset_index()
# replace NaN with 0
df_coop["protein1"].fillna(0,inplace=True)
# rename
df_coop.rename(columns={"protein1":"num_partners","cooperativity":"partner_cooperativity"},inplace=True)


df_tf=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_merged.csv")
# merge df_coop with df_tf, on protein2
df_coop=pd.merge(df_coop,df_tf,left_on="protein2",right_on="protein2",how="inner")
df_coop["protein2_cooperativity"]=pd.cut(df_coop["cooperativity_index"],bins=[-1,0.3,0.7,2],labels=["redundant","other","codependent"])


# plot distribution of num_partners, color by partner_cooperativity
sns.boxplot(x="protein2_cooperativity",
            y="num_partners",
            hue="partner_cooperativity",
            data=df_coop)
plt.savefig("partner_specificity.png")
plt.close()
