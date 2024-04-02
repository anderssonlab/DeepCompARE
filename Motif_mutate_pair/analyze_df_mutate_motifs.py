import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

df=pd.read_csv('df_mutate_pair_promoters_broad_common_remap_Hep-G2,track0.csv',header=None)
df.shape


# group by tf1, calculate average additivity and AND_relation
# reveal general cooperitivity of TFs
df_grouped=df.groupby(["tf1"]).agg({"Additivity":np.mean,"AND_relation":np.mean})


# plot the results
texts = []
plt.figure(figsize=(10, 8))  
scatter_plot = sns.scatterplot(data=df_grouped, x='Additivity', y='AND_relation')
for line in range(0, df_grouped.shape[0]):
    text=scatter_plot.text(df_grouped.Additivity[line], df_grouped.AND_relation[line], 
                            df_grouped.index[line], horizontalalignment='left', 
                            size='small', color='black')
    texts.append(text)
adjust_text(texts)
plt.title("AVG Additivity vs AVG AND_relation")
plt.savefig("additivity_AND_relation.pdf",dpi=300)
plt.close()



# group by tf1 and tf2, calculate average additivity and AND_relation
# then do hierarchical clustering on tfs
# plot heatmap

df_grouped=df.groupby(["tf1","tf2"]).agg({"Additivity":np.mean,"AND_relation":np.mean})
df_grouped=df_grouped.reset_index()
df_grouped=df_grouped.pivot(index="tf1",columns="tf2",values="Additivity")
df_grouped=df_grouped.fillna(0)
plt.figure(figsize=(20,18))
sns.clustermap(df_grouped)
plt.title("AVG Additivity")
plt.savefig("additivity_heatmap.pdf",dpi=300)
plt.close()

df_grouped=df.groupby(["tf1","tf2"]).agg({"Additivity":np.mean,"AND_relation":np.mean})
df_grouped=df_grouped.reset_index()
df_grouped=df_grouped.pivot(index="tf1",columns="tf2",values="AND_relation")
df_grouped=df_grouped.fillna(0)
plt.figure(figsize=(20,18))
sns.clustermap(df_grouped)
plt.title("AVG AND_relation")
plt.savefig("AND_relation_heatmap.pdf",dpi=300)
plt.close()

