import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_cooperativity_index_merged.csv")
df_annot=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors]/DatabaseExtract_v_1.01.csv",index_col=0)
df=pd.merge(df,df_annot,left_on="protein2",right_on="HGNC symbol",how="inner")



df=df[['protein2','cooperativity_index','DBD']].copy()
# order df by cooperativity index
df[df.cooperativity_index>0.7].DBD.value_counts()
df[df.cooperativity_index<0.3].DBD.value_counts()

df["redundancy"]=(df.cooperativity_index<0.3).astype(int)
df["codependency"]=(df.cooperativity_index>0.7).astype(int)
df["C2H2"]=df.DBD.str.contains("C2H2").astype(int)
df['bHLH']=df.DBD.str.contains("bHLH").astype(int)
df["Ets"]=df.DBD.str.contains("Ets").astype(int)
df["Nuclear receptor"]=df.DBD.str.contains("Nuclear receptor").astype(int)
df["bZIP"]=df.DBD.str.contains("bZIP").astype(int)
df["Forkhead"]=df.DBD.str.contains("Forkhead").astype(int)
df["Homeodomain"]=df.DBD.str.contains("Homeodomain").astype(int)

# chi2 test
chi2_contingency(pd.crosstab(df.C2H2,df.redundancy))
chi2_contingency(pd.crosstab(df.Ets,df.redundancy))

chi2_contingency(pd.crosstab(df["Nuclear receptor"],df.codependency))
chi2_contingency(pd.crosstab(df.bZIP,df.codependency))
chi2_contingency(pd.crosstab(df.Forkhead,df.codependency))
chi2_contingency(pd.crosstab(df.Homeodomain,df.codependency))


# select Ets, C2H2, "Nuclear receptor", Forkhead, Homeodomain
df=df[df["DBD"].isin(["Ets","C2H2 ZF", "Nuclear receptor","Forkhead","Homeodomain"])].copy()


sns.boxplot(data=df,x="DBD",y="cooperativity_index")
# rotate x-axis labels
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("dbd_vs_ci.pdf")
plt.close()


