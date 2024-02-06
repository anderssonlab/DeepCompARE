"""
read in res.csv
calculate correlation between AF (column zero) and prediction (column 1-16)
calculate correlation between AF (column zero) and absolute value of prediction (column 1-16)
plot the correlation using bar plot
"""

import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Corr_with_maf/res.csv",header=None)

corrs=[]
pvals=[]
tracks=[]
condition_abs=[]
for i in range(1,17):
    logger.info(f"processing track {i}")
    corr,pval=pearsonr(df.iloc[:,0],df.iloc[:,i])
    corrs.append(corr)
    pvals.append(pval)
    tracks.append(i)
    condition_abs.append(False)
    
    corr,pval=pearsonr(df.iloc[:,0],abs(df.iloc[:,i]))
    corrs.append(corr)
    pvals.append(pval)
    tracks.append(i)
    condition_abs.append(True)
    

df_corr=pd.DataFrame({"tracks":tracks,"corrs":corrs,"pvals":pvals,"condition_abs":condition_abs})
df_corr["tracks"]=df_corr["tracks"].astype(str)
df.to_csv("/isdata/alab/people/pcr980/DeepCompare/Corr_with_maf/corrs.csv",index=False,header=False)

sns.barplot(x="tracks",y="corrs",data=df_corr,hue="condition_abs")
plt.savefig("/isdata/alab/people/pcr980/DeepCompare/Corr_with_maf/corr.png")