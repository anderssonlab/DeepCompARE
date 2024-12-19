import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
# Do TF pair coming from same family have lower cooperativity ratio?


for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
    df=df[df["c_sum"]>5].reset_index(drop=True)
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
    
    # mannwhitneyu test
    print(mannwhitneyu(df[df["same_family"]==True]["cooperativity_index"],df[df["same_family"]==False]["cooperativity_index"]))
    
    # sns kde plot the distribution of cooperativity_index split by family with family size
    sns.kdeplot(
        data=df,
        x='cooperativity_index',
        hue='same_family',
        cut=0,
        common_norm=False,
        fill=True
    )

    plt.savefig(f"same_family_{cell_line}.pdf")
    plt.close()
