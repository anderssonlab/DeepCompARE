import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_tf_hepg2=pd.read_csv('tf_cooperativity_index_hepg2.csv').loc[:,['protein2','cooperativity_index']].copy()
df_tf_k562=pd.read_csv('tf_cooperativity_index_k562.csv').loc[:,['protein2','cooperativity_index']].copy()

# merge by protein2 and rename columns
df_tf_hepg2.columns=['protein2','cooperativity_index_hepg2']
df_tf_k562.columns=['protein2','cooperativity_index_k562']
df_tf=pd.merge(df_tf_hepg2,df_tf_k562,on='protein2',how='inner')

# pearson correlation
df_tf.corr()
# scatter plot with annotation
sns.scatterplot(x='cooperativity_index_hepg2',y='cooperativity_index_k562',data=df_tf)
for i in range(df_tf.shape[0]):
    if abs(df_tf.iloc[i,1]-df_tf.iloc[i,2])>0.3:
        plt.text(df_tf.iloc[i,1],df_tf.iloc[i,2],df_tf.iloc[i,0],fontsize=6)

plt.savefig("Plots_cell_type_specificity/scatter_ci_hepg2_vs_k562.pdf")
plt.close()







df_tf_promoter=pd.read_csv('tf_cooperativity_index_promoter.csv').loc[:,['protein2','cooperativity_index']].copy()
df_tf_enhancer=pd.read_csv('tf_cooperativity_index_enhancer.csv').loc[:,['protein2','cooperativity_index']].copy()
# merge by protein2 and rename columns
df_tf_promoter.columns=['protein2','cooperativity_index_promoter']
df_tf_enhancer.columns=['protein2','cooperativity_index_enhancer']
df_tf=pd.merge(df_tf_promoter,df_tf_enhancer,on='protein2',how='inner')

# pearson correlation
df_tf.corr()
# scatter plot with annotation
sns.scatterplot(x='cooperativity_index_promoter',y='cooperativity_index_enhancer',data=df_tf)
for i in range(df_tf.shape[0]):
    if abs(df_tf.iloc[i,1]-df_tf.iloc[i,2])>0.4:
        plt.text(df_tf.iloc[i,1],df_tf.iloc[i,2],df_tf.iloc[i,0],fontsize=6)

plt.savefig("Plots_cell_type_specificity/scatter_ci_promoter_vs_enhancer.pdf")
plt.close()
