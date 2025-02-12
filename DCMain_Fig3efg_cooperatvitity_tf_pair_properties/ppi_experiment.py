import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

for bait in ["BACH1","RFX5","IKZF1","MAFG","RREB1"]:
    df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_k562_dhs.csv")
    # select only c_sum > 1
    df_coop = df_coop[df_coop["c_sum"] > 1].reset_index(drop=True)
    df_coop["c_redundancy"] = df_coop["c_redundancy"].abs()
    df_ppi = pd.read_csv(f"Pd1_Petra_data/2025-01-09_s10_PublishedPPIandProtComplexes_0.5_0.1_k562.txt", sep='\t')
    df_ppi.rename(columns={'SWI/SNF': 'SWI_SNF'}, inplace=True)
    df_coop = pd.merge(df_coop, df_ppi, on=["protein1", "protein2"], how="left")
    # change  Reported_PPI to 0,1, 0 for no, 1 for yes
    df_coop["Reported_PPI"] = df_coop["Reported_PPI"].fillna(0)
    df_coop["Reported_PPI"] = df_coop["Reported_PPI"].map({"No": 0, "Yes": 1})
    df_coop = df_coop.dropna(subset=["Reported_PPI"]).reset_index(drop=True)
    df_coop_bait=df_coop[df_coop["protein2"]==bait].reset_index(drop=True)
    df_exp=pd.read_csv(f"Pd2_PPI_experiment/250107.{bait}.K562.Taplin.GenoppiStats.txt", sep='\t')
    # which genes in df_exp are in df_coop_bait?
    candidates=df_coop_bait["protein1"].tolist()
    # split dimers in candidates
    candidates = [x.split("::") for x in candidates]
    # flatten the list
    candidates = [item for sublist in candidates for item in sublist]


    tfs_found=set(df_exp["gene"].values).intersection(candidates)
    df_coop_bait_found = df_coop_bait[df_coop_bait["protein1"].isin(tfs_found)]





    # Plot and annotate
    fig, axes = plt.subplots(3, 2, figsize=(6, 8), sharex=True)
    # Define the columns to plot, including 'distance'
    columns_to_plot = ['cooperativity_index', 'c_sum', 'c_redundancy', 'c_codependency', 'distance']
    titles = ['Cooperativity Index', 'C Sum', 'C Redundancy', 'C Codependency', 'Distance']
    # Adjust axes for a 3x2 grid
    axes = axes.flatten()
    # Plot each column in a subplot
    for ax, column, title in zip(axes, columns_to_plot, titles):
        # decrease the size of the dots
        sns.swarmplot(ax=ax, x="Reported_PPI", y=column, data=df_coop_bait, size=3)
        for _, row in df_coop_bait_found.iterrows():
            ax.text(
                x=row["Reported_PPI"],
                y=row[column],
                s=row["protein1"],
                color="red",
                ha="center"
            )
        ax.set_title(title)
        ax.set_xlabel("Reported PPI")
        ax.set_ylabel(column)

    #
    # Remove empty subplot if the number of plots is less than the grid size
    if len(columns_to_plot) < len(axes):
        for ax in axes[len(columns_to_plot):]:
            ax.remove()


    plt.tight_layout()
    plt.savefig(f"{bait}.png")
    plt.close()
    
    
    
# coop index v.s. log fold change