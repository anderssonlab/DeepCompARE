import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

from scipy.stats import pearsonr
import numpy as np

# TODO: add Pol II ChIP-Seq

re="proximal"




# read histone marks
df_histone=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd10_chromatin_profile/{re}_k562.csv")
# add region in format chr:start-end
df_histone["region"]=df_histone["chrom"].astype(str)+":"+df_histone["start"].astype(str)+"-"+df_histone["end"].astype(str)
# remove chrom,start,end
df_histone=df_histone.drop(columns=["chrom","start","end"])


# read mediator interactors
df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)
# output the "gene" in text format
# pd.Series(np.sort(df_med["gene"].unique())).to_csv("mediator_interactors.txt",index=False,header=False)



# read chip
df_chip=pd.read_csv(f"chip_{re}.csv")
df_chip["region"]=df_chip["chrom"].astype(str)+":"+df_chip["start"].astype(str)+"-"+df_chip["end"].astype(str)
df_chip=df_chip.drop(columns=["chrom","start","end"])
# merge df_chip with df_histone
df=pd.merge(df_chip,df_histone,on="region",how="inner")



for mediator in df_med["bait"].unique().tolist()+["all"]:
    if mediator=="all":
        mediator_interactors=df_med["gene"].unique().tolist()
    else:
        mediator_interactors=df_med[df_med["bait"]==mediator]["gene"].unique().tolist()
    mediator_interactor_cols=[f"chip_{tf}" for tf in mediator_interactors]
    # select columns in df
    mediator_interactor_cols=[col for col in mediator_interactor_cols if col in df.columns]
    df[f"sum_chip_{mediator}_interactors (n={len(mediator_interactor_cols)})"]=df[mediator_interactor_cols].sum(axis=1)

# remove columns with all 0
df=df.loc[:,(df!=0).any(axis=0)]



def set_insigificant_to_zero(df,df_corr):
        # if p_value>0.05, set to 0
    for rowname in df_corr.index.tolist():
        for colname in df_corr.columns.tolist():
            r,p=pearsonr(df[rowname],df[colname])
            if p>0.05:
                logger.info(f"set {rowname} and {colname} to 0")
                df_corr.loc[rowname,colname]=0
    return df_corr



correlations=df.corr()
corr_histone_sum_chip=correlations.loc[correlations.index.str.startswith("sum_chip"),correlations.columns.str.startswith("log1p")]
corr_histone_sum_chip=set_insigificant_to_zero(df,corr_histone_sum_chip)

# clustermap
sns.clustermap(corr_histone_sum_chip,annot=True,figsize=(8,6),fmt=".2f",cmap="coolwarm",vmin=-0.6,vmax=0.6)
plt.title(f"{re}")
plt.tight_layout()
plt.savefig(f"chip_corr_histone_sum_chip_{re}.pdf")
plt.close()


corr_histone_tf_chip=correlations.loc[correlations.index.str.startswith("log1p"),correlations.columns.str.startswith("chip_")]
corr_histone_tf_chip=set_insigificant_to_zero(df,corr_histone_tf_chip)
# clustermap
sns.clustermap(corr_histone_tf_chip,annot=True,figsize=(18,8),fmt=".2f",cmap="coolwarm",vmin=-0.6,vmax=0.6)
plt.title(f"{re}")
plt.tight_layout()
plt.savefig(f"chip_corr_histone_tf_chip_{re}.pdf")
plt.close()

