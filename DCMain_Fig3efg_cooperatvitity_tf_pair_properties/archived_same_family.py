import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
# Do TF pair coming from same family have lower cooperativity ratio?


for cell_line in ["hepg2","k562","merged"]:
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

    # retain the "protein2" with at least 40 pair with same_family==True
    df_filtered=pd.DataFrame()
    for family,df_sub in df.groupby("family_2"):
        if df_sub["same_family"].sum()>=40:
            print(family)
            print(df_sub["same_family"].value_counts())
            df_filtered=pd.concat([df_filtered,df_sub],axis=0)

    # sns violin plot the distributon of cooperativity_index split by family
    # Calculate the size of each family in family_2
    family_sizes = df_filtered["family_2"].value_counts().reset_index()
    family_sizes.columns = ["family_2", "family_size"]

    # Merge the family sizes back into the filtered DataFrame
    df_filtered = pd.merge(df_filtered, family_sizes, on="family_2", how="left")

    # sns violin plot the distribution of cooperativity_index split by family with family size
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        data=df_filtered,
        x='family_2',
        y='cooperativity_index',
        hue="same_family",
        split=True,
        cut=0,
        quantiles=[0.25, 0.5, 0.75],
        inner=None
    )
    plt.xticks(rotation=90)

    # Annotate the x-axis labels with family size
    xticks_labels = [f"{family}\n(n={size})" for family, size in zip(family_sizes["family_2"], family_sizes["family_size"])]
    plt.gca().set_xticklabels(xticks_labels)

    # Make bottom margin larger
    plt.gcf().subplots_adjust(bottom=0.6)
    plt.savefig(f"same_family_{cell_line}.pdf")
    plt.close()
