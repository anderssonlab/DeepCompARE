import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

#---------------------------
# 1. Effect of nonlinearity weight: compare cooperativity_index and cooperativity_fraction
#---------------------------
fig,axs=plt.subplots(2,3,figsize=(12,8))
for i,suffix in enumerate(["hepg2","k562","merged"]):
    for j,data_type in enumerate(["tf","tf_pair"]):
        ax=axs[j,i]
        df=pd.read_csv(f'{data_type}_cooperativity_index_{suffix}.csv')
        df=df[df['c_sum']>1].reset_index(drop=True)
        ax.scatter(df['cooperativity_index'],df['cooperativity_fraction'],s=1)
        r=pearsonr(df['cooperativity_index'],df['cooperativity_fraction'])[0]
        ax.text(0.5,0.1,f'Pearson r={r:.2f}',transform=ax.transAxes)
        ax.set_title(f'{data_type} level, {suffix}')
        if j==1:
            ax.set_xlabel('cooperativity_index (weighted)')
        if i==0:
            ax.set_ylabel('cooperativity_fraction (unweighted)' )

plt.savefig('compare_cooperativity_index_vs_fraction.pdf')
plt.close()



#---------------------------
# 2. cell type effect
#---------------------------
fig,axs=plt.subplots(1,2,figsize=(12,4))
for i,data_type in enumerate(["tf","tf_pair"]):
    ax=axs[i]
    df_hepg2=pd.read_csv(f'{data_type}_cooperativity_index_hepg2.csv')
    df_hepg2=df_hepg2[df_hepg2['c_sum']>1].reset_index(drop=True)
    df_k562=pd.read_csv(f'{data_type}_cooperativity_index_k562.csv')
    df_k562=df_k562[df_k562['c_sum']>1].reset_index(drop=True)
    # columns to merge is all columns starting with protein
    col_merge=[col for col in df_hepg2.columns if col.startswith('protein')]
    df_merged=pd.merge(df_hepg2,df_k562,on=col_merge,suffixes=('_hepg2','_k562'))
    ax.scatter(df_merged['cooperativity_index_hepg2'],df_merged['cooperativity_index_k562'],s=1)
    r=pearsonr(df_merged['cooperativity_index_hepg2'],df_merged['cooperativity_index_k562'])[0]
    ax.text(0.5,0.1,f'Pearson r={r:.2f}',transform=ax.transAxes)
    ax.set_title(f'{data_type} level')
    ax.set_xlabel('cooperativity_index (HepG2)')
    if i==0:
        ax.set_ylabel('cooperativity_index (K562)' )

plt.savefig('compare_hepg2_vs_k562.pdf')
plt.close()


#---------------------------
# 3. regulatory element effect
#---------------------------
fig,axs=plt.subplots(1,2,figsize=(12,4))
for i,data_type in enumerate(["tf","tf_pair"]):
    ax=axs[i]
    df_promoter=pd.read_csv(f'{data_type}_cooperativity_index_promoter.csv')
    df_promoter=df_promoter[df_promoter['c_sum']>1].reset_index(drop=True)
    df_enhancer=pd.read_csv(f'{data_type}_cooperativity_index_enhancer.csv')
    df_enhancer=df_enhancer[df_enhancer['c_sum']>1].reset_index(drop=True)
    # columns to merge is all columns starting with protein
    col_merge=[col for col in df_promoter.columns if col.startswith('protein')]
    df_merged=pd.merge(df_promoter,df_enhancer,on=col_merge,suffixes=('_promoter','_enhancer'))
    ax.scatter(df_merged['cooperativity_index_promoter'],df_merged['cooperativity_index_enhancer'],s=1)
    r=pearsonr(df_merged['cooperativity_index_promoter'],df_merged['cooperativity_index_enhancer'])[0]
    ax.text(0.5,0.1,f'Pearson r={r:.2f}',transform=ax.transAxes)
    ax.set_title(f'{data_type} level')
    ax.set_xlabel('cooperativity_index (promoter)')
    if i==0:
        ax.set_ylabel('cooperativity_index (enhancer)' )

plt.savefig('compare_promoter_vs_enhancer.pdf')
plt.close()


