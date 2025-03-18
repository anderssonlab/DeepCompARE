import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_hepg2=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_hepg2_pe.csv")
df_hepg2["cell"]="hepg2"
df_k562=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_k562_pe.csv")
df_k562["cell"]="k562"
df=pd.concat([df_hepg2,df_k562],axis=0)


df_annot=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv",index_col=0)
df=pd.merge(df,df_annot,left_on="protein2",right_on="HGNC symbol",how="inner")
df=df[['protein2','cooperativity_index','DBD','cell']].copy()
df_counts=df.DBD.value_counts()
dbds_retain=df_counts[df_counts>10].index
# filter df by DBD, retain only DBDs with at least 5 entry
df=df[df.DBD.isin(dbds_retain)].reset_index(drop=True)

sns.boxplot(data=df,x="DBD",y="cooperativity_index",hue="cell",palette="Set2")
# rotate x-axis labels
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("dbd_vs_ci.pdf")
plt.close()





