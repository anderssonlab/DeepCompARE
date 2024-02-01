import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

df=pd.read_csv('df_mutate_motif_cage_k562.csv',header=None,nrows=10000000)
df.shape
df.columns=["idx1","idx2","region_idx",
            "chromosome","start1","end1",
            "tf1","start_rel1","end_rel1",
            "chromosome2","start2","end2",
            "tf2","start_rel2","end_rel2",
            "pred_orig","pred_mut1","pred_mut2","pred_mut_both",
            "ism_score_mut1","ism_score_mut2","ism_score_mut_both",
            "AND_relation","Additivity"]
# check the following columns are almost equal
np.allclose(df.ism_score_mut_both,df.pred_orig-df.pred_mut_both)
np.allclose(df.ism_score_mut1,df.pred_orig-df.pred_mut1)
np.allclose(df.ism_score_mut2,df.pred_orig-df.pred_mut2)
np.allclose(df.Additivity,df.ism_score_mut1+df.ism_score_mut2-df.ism_score_mut_both)
np.allclose(df.AND_relation,(df.ism_score_mut1+df.ism_score_mut2)/2-df.ism_score_mut_both)

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

