import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt





import matplotlib
matplotlib.rcParams['pdf.fonttype']=42





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

# sort each group by cooperativity_index,
order=df.groupby("DBD")["cooperativity_index"].mean().sort_values(ascending=False).index
df["DBD"] = pd.Categorical(df["DBD"], categories=order, ordered=True)

plt.figure(figsize=(3,4.5))
# thin frame, linewidth 0.5
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
# small fliers
sns.boxplot(data=df,y="DBD",x="cooperativity_index",hue="cell",palette="Set2",linewidth=0.5,fliersize=1)
# rotate x-axis labels
plt.ylabel("DNA-binding domain (DBD)", fontsize=7)
plt.xlabel("TF synergy score", fontsize=7)
plt.xticks(rotation=90, fontsize=5)
plt.yticks(fontsize=5)
plt.legend(fontsize=5, title_fontsize=5)
plt.tight_layout()
plt.savefig("dbd_vs_ci.pdf")
plt.close()





