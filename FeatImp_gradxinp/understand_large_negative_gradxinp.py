import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns   
from scipy.stats import pearsonr

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/motif_df.csv")
# remove rows with nan
df=df[~np.isnan(df["max_gradxinp"])]

corr,pval=pearsonr(df["max_gradxinp"],df["min_gradxinp"]) # -0.82,0.0
corr,pval=pearsonr(df["max_gradxinp"],df["min_gradxinp"].abs()) # 0.82, 0
corr,pval=pearsonr(df["max_gradxinp"],df["mean_gradxinp"]) # 0.55, 0
corr,pval=pearsonr(df["max_gradxinp"],df["mean_abs_gradxinp"]) # 0.94, 0
corr,pval=pearsonr(df["max_gradxinp"],df["max_abs_gradxinp"]) # 0.97,0
corr,pval=pearsonr(df["max_abs_gradxinp"],df["mean_abs_gradxinp"]) # 0.97, 0

# randomly sample 1000 rows from df
df_sub=df.sample(n=1000,random_state=1)

# scatter plot max_gradxinp vs min_gradxinp

sns.scatterplot(x="max_gradxinp",y="min_gradxinp",data=df_sub)
plt.title("Randomly sampled 1000")
# add correlation and pvalue
plt.text(0.5,0.5,"corr: "+str(round(corr,2))+"\npval: "+str(round(pval,2)))
plt.legend()
plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Putative_binding_vs_actual_binding/Plots_understand_gradxinp/max_gradxinp_vs_min_gradxinp.png",dpi=300)
plt.close()
