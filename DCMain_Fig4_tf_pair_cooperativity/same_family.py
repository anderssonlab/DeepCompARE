import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
# Do TF pair coming from same family have lower cooperativity ratio?

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}_pe.csv")
    df=assign_cooperativity(df,1,0.9,0.3,0.7)
    df=df[df["cooperativity"]!="Linear"].reset_index(drop=True)
    df_family=pd.read_csv("/isdata/alab/people/pcr980/Resource/JASPAR2022_tracks/JASPAR2024_CORE_tf_family.csv")
    df_family["ID"]=df_family["ID"].str.upper()

    # are all df["protein2"] in df_family["ID"]?
    df=pd.merge(df,df_family,left_on="protein1",right_on="ID",how="inner")
    df.drop(columns=['AC', 'ID'],inplace=True)
    df.rename(columns={"tf_family":"family_1"},inplace=True)
    df=pd.merge(df,df_family,left_on="protein2",right_on="ID",how="inner")
    df.drop(columns=['AC', 'ID'],inplace=True)
    df.rename(columns={"tf_family":"family_2"},inplace=True)
    df["same_family"]= (df["family_1"]==df["family_2"])
    
    
    # sns kde plot the distribution of cooperativity_index split by family with family size
    plt.figure(figsize=(2.3, 2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.kdeplot(
        data=df,
        x='cooperativity_index',
        hue='same_family',
        cut=0,
        common_norm=False,
        fill=True,
        linewidth=0.5,
        color="coolwarm",
    )
    # write the p value to the plot
    stat, p= mannwhitneyu(df[df["same_family"]==True]["cooperativity_index"],df[df["same_family"]==False]["cooperativity_index"])
    plt.text(0.1, 0.2, f"p={p:.2e}", fontsize=7, transform=plt.gca().transAxes)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel("Cooperativity index", fontsize=7)
    plt.ylabel("Density", fontsize=7)
    plt.legend(
    title="Same TF family?",
    fontsize=5,
    title_fontsize=5,
    loc='upper right',
    labels=["True", "False"]  # Explicitly set legend labels
)
    plt.tight_layout()
    plt.savefig(f"same_family_{cell_line}.pdf")
    plt.close()