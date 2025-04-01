import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt




# read gini index
df_gini=pd.read_csv('/isdata/alab/people/pcr980/Resource/gtex.dispersionEstimates.tab',sep='\t')

# read mediator interactors
df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
df_med=df_med[df_med["significant"]].reset_index(drop=True)

# read cofactors
df_cof=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/CofactorComplexes.txt",sep="\t")

# read human transcription factors
def get_htfs_list():
    # read human transcription factors
    df_htfs=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/Human_transcription_factors/DatabaseExtract_v_1.01.csv",index_col=0)
    # subset for "Is TF?"=="Yes" or "HGNC symbol"=="SMAD2"
    df_htfs=df_htfs[(df_htfs["Is TF?"]=="Yes")|(df_htfs["HGNC symbol"]=="SMAD2")].reset_index(drop=True)
    return df_htfs["HGNC symbol"].tolist()
    

tf_list=get_htfs_list()




df_gini["is_med_interactor"]=df_gini["symbol"].isin(df_med["gene"].unique())
df_gini["is_cofactor"]=df_gini["symbol"].isin(df_cof["Gene"].unique())
df_gini["is_tf"]=df_gini["symbol"].isin(tf_list)
df_gini["non_med_interactor_tf"]=df_gini["is_tf"] & ~df_gini["is_med_interactor"]





plt.figure(figsize=(2.5,2))
# thin frame
plt.gca().spines['top'].set_linewidth(0.5)
plt.gca().spines['right'].set_linewidth(0.5)
plt.gca().spines['bottom'].set_linewidth(0.5)
plt.gca().spines['left'].set_linewidth(0.5)
sns.kdeplot(df_gini[df_gini["is_med_interactor"]]["gini"],label="Mediator-Interacting TFs",linewidth=1)
sns.kdeplot(df_gini[df_gini["is_cofactor"]]["gini"],label="Cofactor",linewidth=1)
sns.kdeplot(df_gini[~df_gini["is_med_interactor"] & ~df_gini["is_cofactor"]]["gini"],label="Others",linewidth=1)
sns.kdeplot(df_gini[df_gini["non_med_interactor_tf"]]["gini"],label="Non-Mediator-Interacting TFs",linewidth=1)
plt.legend(fontsize=5,title_fontsize=5)
plt.xlabel("Gini Index of Gene Expression",fontsize=7)
plt.ylabel("Density",fontsize=7)
plt.xticks(fontsize=5)
plt.yticks(fontsize=5)
plt.tight_layout()
plt.savefig("gini_index.pdf")
plt.close()











#---------------------------------------
# can you predict TF gene specificity given DNA sequence?
#---------------------------------------

df_gini=df_gini[df_gini["is_tf"]].reset_index(drop=True)
