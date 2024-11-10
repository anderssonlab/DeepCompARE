import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, confusion_matrix
from scipy.stats import ks_2samp,pearsonr


#--------------------
# Helper functions
#--------------------

# def add_expr_info(df,col):
#     df_expr=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/RNA_expression/Raw_data/ENCFF928NYA_RNASeq_K562.tsv",sep='\t')
#     # remove .xxx from gene_id
#     df_expr["gene_id"]=df_expr["gene_id"].str.split(".").str[0]
#     df_name_mapper=pd.read_csv("/isdata/alab/people/pcr980/Resource/Human_transcription_factors]/DatabaseExtract_v_1.01.csv")
#     df_name_mapper=df_name_mapper[['Ensembl ID', 'HGNC symbol']].copy()
#     # merge with name mapper
#     df_expr=pd.merge(df_expr,df_name_mapper,left_on='gene_id',right_on='Ensembl ID',how='inner')
#     df_expr=df_expr[['HGNC symbol', 'TPM']].copy()
#     # merge with df
#     df=pd.merge(df,df_expr,left_on=col,right_on='HGNC symbol',how='inner')
#     df.drop(columns=['HGNC symbol'],inplace=True)
#     return df


def plot_col_by_hue(df,col,hue):
        df_sub=df[[hue,col]].copy()
        df_sub[hue]=df_sub[hue].astype('category')
        # remove NaN
        df_sub=df_sub.dropna()
        sns.kdeplot(data=df_sub,x=col,hue=hue,common_norm=False)
        # mann whitney u test, allow NaN
        # get unique values in hue
        unique_values=df_sub[hue].unique()
        e=0
        for i in range(len(unique_values)):
            for j in range(i+1,len(unique_values)):
                u_stat,p=mannwhitneyu(df_sub[df_sub[hue]==unique_values[i]][col],df_sub[df_sub[hue]==unique_values[j]][col],alternative='two-sided')
                plt.text(df_sub[col].max(),0.5+e*0.1,f"{unique_values[i]} vs {unique_values[j]} p={p:.2e}")
                e=e+1
        plt.title(hue)
        plt.savefig(f'Plots/{col}_vs_{hue}.pdf')
        plt.close()



def signed_ks_test(df,col):
    dstat, p = ks_2samp(df[df["Reported_PPI_lowerTh"]=="Yes"][col].abs(),df[df["Reported_PPI_lowerTh"]=="No"][col].abs())
    if df[df["Reported_PPI_lowerTh"]=="Yes"][col].abs().median()<df[df["Reported_PPI_lowerTh"]=="No"][col].abs().median():
        dstat=-dstat
    return dstat, p



#---------------------------------
# Read data and merge
#---------------------------------
df=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_cooperativity/tf_pair_cooperativity_index_k562.csv")
df["count_sum"]=df["count_codependency"]+df["count_redundancy"]
df=df[df["count_sum"]>5].reset_index(drop=True)

df_ppi=pd.read_csv("2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt",sep='\t')

df_ppi=df_ppi[['protein1', 'protein2', 'Reported_PPI_lowerTh']].copy()
df_ppi2=df_ppi.copy()
df_ppi2.rename(columns={'protein1':'protein2','protein2':'protein1'},inplace=True)
df_ppi=pd.concat([df_ppi,df_ppi2],axis=0).reset_index(drop=True)

df=pd.merge(df,df_ppi,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')



# # add expr_info if necessary
# df=add_expr_info(df,'protein1')
# df.rename(columns={'TPM':'TPM1'},inplace=True)
# df=add_expr_info(df,'protein2')
# df.rename(columns={'TPM':'TPM2'},inplace=True)
# df["min_expr"]=df[['TPM1','TPM2']].min(axis=1)
# df["mean_expr"]=df[['TPM1','TPM2']].mean(axis=1)




#---------------------------------
# Plot distribution split by Reported_PPI_lowerTh
#---------------------------------

for col in ['c_codependency','count_codependency', 'distance', 'distance_iqr','c_sum', 'cooperativity_index', 'cooperativity_fraction', 'count_sum','min_expr', 'mean_expr']:
    hue="Reported_PPI_lowerTh"
    plot_col_by_hue(df,col,hue)



df_ppi=pd.read_csv("2024-10-24_s1_PublishedPPIannotation_withProtComplexes.txt",sep='\t')
# rename 'SWI/SNF_lowerTh' to 'SWI_SNF_lowerTh'
df_ppi.rename(columns={'SWI/SNF_lowerTh':'SWI_SNF_lowerTh',
                       'SWI/SNF':'SWI_SNF',
                       },inplace=True)
for hue in ['Mediator', 'TFIIB', 'TFIID', 'TFIIE', 'TFIIF', 'TFIIH', 'SWI_SNF', 'POLII']:
    plot_col_by_hue(df_ppi,'coopIndex_avg',hue)




#-----------------------------
# analyze chip-seq colocolization info
#-----------------------------
df_chip_coloc=pd.read_csv("/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_K562_summary.csv")
# rename tf1 to protein1, tf2 to protein2
df_chip_coloc.rename(columns={'tf1':'protein1','tf2':'protein2'},inplace=True)
df=pd.merge(df,df_chip_coloc,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')
df.rename(columns={'mean_abs':'mean_chip_peak_distance',
                   'tf_pair_count':'chip_overlap_count',
                   },inplace=True)

plot_col_by_hue(df,'chip_overlap_count','Reported_PPI_lowerTh')
plot_col_by_hue(df,'mean_chip_peak_distance','Reported_PPI_lowerTh')











#-----------------------------
# random forest
#-----------------------------
def train_rf(features):
    df_sub = df[features].copy()
    # Drop NaN values and reset the index
    df_sub = df_sub.dropna().reset_index(drop=True)
    X = df_sub.drop(columns='Reported_PPI').values
    Y = df_sub["Reported_PPI"].values
    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    # Train the Logistic Regression model
    clf = RandomForestClassifier(n_estimators=10)
    clf.fit(X_train, y_train)
    # Predict with the model
    y_pred = clf.predict(X_test)
    # Check the accuracy
    f1 = f1_score(y_test, y_pred)
    conf_matrix = confusion_matrix(y_test, y_pred)
    return f1, conf_matrix



train_rf(['Reported_PPI', 'cooperativity_index'])
train_rf(['Reported_PPI', 'c_codependency'])
train_rf(['Reported_PPI', 'chip_overlap_count'])
train_rf(['Reported_PPI', 'mean_chip_peak_distance'])
train_rf(['Reported_PPI', 'c_codependency','c_redundancy'])
train_rf(['Reported_PPI', 'cooperativity_index','mean_chip_peak_distance'])
train_rf(['Reported_PPI', 'cooperativity_index','distance_codependency','distance_redundancy','count_codependency','count_redundancy'])





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


for tf in df['protein2'].unique():
    df_sub=df[df['protein2']==tf].copy()
    if df_sub["Reported_PPI_lowerTh"].nunique()<2:
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
df_res[(df_res['ci_stat']<0) & (df_res['ci_p']<0.05)]
df_outliers=df_res[(df_res['ci_stat']<0) & (df_res['ci_p']<0.05)].reset_index(drop=True)
df_outliers.to_csv("outliers.csv",index=False)
df_outliers[(df_outliers['c_redundancy_stat']>0) & (df_outliers['c_redundancy_p']<0.05)].reset_index(drop=True)
df_outliers[(df_outliers['c_codependency_stat']<0)].reset_index(drop=True)
df_outliers.tf.values


# for each tf in df_outliers,
# plot the distribution of cooperativity_index, c_codependency, c_redundancy
for tf in df_outliers['tf'].unique():
    pass
    
tf="USF1"
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





df_dispersion=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Cooperativity_effects_tf_level/TFs.dispersionEstimates.k562.tab",sep="\t")
df_dispersion=pd.read_csv("/isdata/alab/people/pcr980/DeepCompare/Cooperativity_effects_tf_level/TFs.dispersionEstimates.hepG2.tab",sep="\t")



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


