import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ks_2samp
from matplotlib.backends.backend_pdf import PdfPages
from loguru import logger



#---------------------------------
# Read and merge data
#---------------------------------
df_coop=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562.csv")
df_coop=df_coop[df_coop['c_sum']>1].reset_index(drop=True)
df_ppi=pd.read_csv("/isdata/alab/people/pcr980/Resource/2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt",sep='\t')
df_ppi.rename(columns={'SWI/SNF_lowerTh':'SWI_SNF_lowerTh',
                       'SWI/SNF':'SWI_SNF',
                       },inplace=True)
# remove columns protPair','coopIndex_avg', 'sum_ci_avg'
df_ppi.drop(columns=['protPair','coopIndex_avg', 'sum_ci_avg'],inplace=True)
df_ppi2=df_ppi.copy()
df_ppi2.rename(columns={'protein1':'protein2','protein2':'protein1'},inplace=True)
df_ppi=pd.concat([df_ppi,df_ppi2],axis=0).reset_index(drop=True)

df=pd.merge(df_coop,df_ppi,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')


# Create a single PDF file for all plots
with PdfPages('ci_combined.pdf') as pdf:
    for tf in df.protein2.unique():
        df_sub = df[df['protein2'] == tf].reset_index(drop=True)
        counts = df_sub['Reported_PPI_lowerTh'].value_counts()
        count_yes = counts.get('Yes', 0)
        count_no = counts.get('No', 0)
        # Perform Mann-Whitney U test
        group_yes = df_sub[df_sub['Reported_PPI_lowerTh'] == 'Yes']['cooperativity_index']
        group_no = df_sub[df_sub['Reported_PPI_lowerTh'] == 'No']['cooperativity_index']
        if len(group_yes) > 0 and len(group_no) > 0:
            stat, p_value = mannwhitneyu(group_yes, group_no, alternative='two-sided')
        else:
            p_value = float('nan')  # Handle cases with insufficient data
        # Generate the violin plot
        sns.violinplot(
            data=df_sub, 
            x='Reported_PPI_lowerTh', 
            y='cooperativity_index', 
            cut=0,
            order=['No', 'Yes']  # Ensure the order is "No", "Yes"
        )
        # Annotate the plot with the p-value
        plt.title(f"{tf}\nYes: {count_yes}, No: {count_no}")
        plt.text(
            x=0.5, y=max(df_sub['cooperativity_index']) * 0.9,  # Position the text
            s=f'p = {p_value:.3e}', 
            ha='center', fontsize=10
        )
        # Add the current plot to the PDF
        pdf.savefig()  # Save the current figure to the PDF
        plt.close()  # Close the current plot

#-----------------------------
# ks_test on c_sum and cooperativity_index
#-----------------------------

# how many TF have significantly larger cooperativity index and larger c_sum in PPI vs non-PPI 
# using ks-test


tf_list=[]
ci_stat_list=[]
ci_p_list=[]
c_codependency_stat_list=[]
c_codependency_p_list=[]
c_redudancy_stat_list=[]
c_redudancy_p_list=[]

def signed_ks_test(df,col):
    yes_values = df[df['Reported_PPI_lowerTh'] == 'Yes'][col]
    no_values = df[df['Reported_PPI_lowerTh'] == 'No'][col]
    stat, p_value = ks_2samp(yes_values, no_values)
    if yes_values.median() < no_values.median():
        stat = -stat
    return stat, p_value


for tf in df['protein2'].unique():
    df_sub=df[df['protein2']==tf].copy()
    if df_sub["Reported_PPI"].nunique()<2:
        continue
    stat, p = signed_ks_test(df_sub,'cooperativity_index')
    tf_list.append(tf)
    ci_stat_list.append(stat)
    ci_p_list.append(p)
    c_codependency_stat, c_codependency_p = signed_ks_test(df_sub,'c_codependency')
    c_codependency_stat_list.append(c_codependency_stat)
    c_codependency_p_list.append(c_codependency_p)
    c_redudancy_stat, c_redudancy_p = signed_ks_test(df_sub,'c_redundancy')
    c_redudancy_stat_list.append(c_redudancy_stat)
    c_redudancy_p_list.append(c_redudancy_p)


df_res=pd.DataFrame({'tf':tf_list,
                     'ci_stat':ci_stat_list,
                     'ci_p':ci_p_list,
                     'c_codependency_stat':c_codependency_stat_list,
                     'c_codependency_p':c_codependency_p_list,
                     "c_redundancy_stat":c_redudancy_stat_list,
                     "c_redundancy_p":c_redudancy_p_list
                     })




# how many TFs have neg ci_stat and are significant?
df_res[(df_res['ci_stat']<0)]
df_res[(df_res['ci_stat']<0) & (df_res['ci_p']<0.05)]
df_outliers=df_res[(df_res['ci_stat']<0) & (df_res['ci_p']<0.05)].reset_index(drop=True)
df_outliers.to_csv("outliers.csv",index=False)
df_outliers[(df_outliers['c_redundancy_stat']>0) & (df_outliers['c_redundancy_p']<0.05)].reset_index(drop=True)
df_outliers[(df_outliers['c_codependency_stat']<0)].reset_index(drop=True)
df_outliers.tf.values



    
"USF1" in df_coop.protein2




df_sub=df[df['protein2']==tf].reset_index(drop=True)
df_sub["c_redundancy"]=df_sub["c_redundancy"].abs()
df_sub["Reported_PPI_lowerTh"]=df_sub["Reported_PPI_lowerTh"].astype('category')
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
sns.kdeplot(data=df_sub,x='cooperativity_index',ax=axs[0],hue='Reported_PPI_lowerTh',common_norm=False)
sns.kdeplot(data=df_sub,x='c_codependency',ax=axs[1],hue='Reported_PPI_lowerTh',common_norm=False)
sns.kdeplot(data=df_sub,x='c_redundancy',ax=axs[2],hue='Reported_PPI_lowerTh',common_norm=False)
plt.suptitle(tf)
# annotate number of samples by Reported_PPI_lowerTh
n_yes=df_sub[df_sub['Reported_PPI_lowerTh']=='Yes'].shape[0]
n_no=df_sub[df_sub['Reported_PPI_lowerTh']=='No'].shape[0]    
plt.suptitle(f"{tf}, n_yes={n_yes}, n_no={n_no}")
plt.savefig(f'k562_kde_{tf}.pdf')
plt.close()




sns.kdeplot(data=df_dispersion,x='gini')
i=0
for tf in df_outliers['tf'].unique():
    df_sub=df_dispersion[df_dispersion['gene']==tf].copy()
    if df_sub.shape[0]==0:
        continue
    plt.axvline(x=df_sub['gini'].values[0], color='r', linestyle='--')
    # add name of the gene, vertically label
    plt.text(df_sub['gini'].values[0],i,tf,rotation=90)
    i=i+0.3

# are there significant difference in gini index between df_outliers['tf'].unique() and the rest?
df_outliers_gini=df_dispersion[df_dispersion['gene'].isin(df_outliers['tf'].unique())].copy()
df_rest_gini=df_dispersion[~df_dispersion['gene'].isin(df_outliers['tf'].unique())].copy()
stat, p = ks_2samp(df_outliers_gini['gini'],df_rest_gini['gini'])
print(df_outliers_gini['gini'].median(),df_rest_gini['gini'].median())
plt.text(df_dispersion['gini'].max(),0.5,f"p={p:.2e}")
plt.title("Gini index for cell type specificity")
plt.savefig('Plots_outliers/gini_outliers.pdf')
plt.close()









df_tf=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_cooperativity_index_hepg2.csv")
df_tf["c_redundancy"]=df_tf["c_redundancy"].abs()




col='c_redundancy'

sns.kdeplot(data=df_tf,x=col)
i=0
for tf in df_outliers['tf'].unique():
    df_sub=df_tf[df_tf['protein2']==tf].copy()
    if df_sub.shape[0]==0:
        continue
    plt.axvline(x=df_sub[col].values[0], color='r', linestyle='--')
    # add name of the gene, vertically label
    plt.text(df_sub[col].values[0],i,tf,rotation=90)
    i=i+0.04

# are there significant difference in col between df_outliers['tf'].unique() and the rest?
df_outliers_col=df_tf[df_tf['protein2'].isin(df_outliers['tf'].unique())].copy()
df_rest_col=df_tf[~df_tf['protein2'].isin(df_outliers['tf'].unique())].copy()
stat, p = ks_2samp(df_outliers_col[col],df_rest_col[col])
print(df_outliers_col[col].median(),df_rest_col[col].median())
plt.text(df_tf[col].max(),0.0001,f"p={p:.2e}")
plt.title(f"TF level {col}")
plt.savefig(f'Plots_outliers/tf_{col}.pdf')
plt.close()


