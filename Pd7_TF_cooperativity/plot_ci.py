import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

for cell_type in ["hepg2","k562"]:
    # 1. histogram of tf_pair ci
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_type}_dhs.csv")
    df=df[df["c_sum"]>1].reset_index(drop=True)
    sns.histplot(df['cooperativity_index'], bins=100)
    plt.xlabel('TF pair cooperativity index')
    plt.ylabel('Frequency')
    plt.title('TF pair cooperativity index distribution')
    plt.savefig(f'plot_hist_tf_pair_ci_{cell_type}_dhs.pdf')
    plt.close()


    # # 2. heatmap of tf_pair ci
    # fig=plt.figure(figsize=(50,50))
    # df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
    # gs = gridspec.GridSpec(1, 2, width_ratios=[0.05, 1])  # Adjust width_ratios to control size of the heatmap and color bar
    # # Create axes for the heatmap and color bar
    # cbar_ax = fig.add_subplot(gs[0, 0])  # Color bar on the left
    # heatmap_ax = fig.add_subplot(gs[0, 1])  # Heatmap on the right
    # sns.heatmap(df, cmap="coolwarm", vmin=0, vmax=1, ax=heatmap_ax, cbar_ax=cbar_ax)
    # plt.margin=0.01
    # plt.savefig(f'plot_heatmap_tf_pair_ci_{cell_type}.pdf')
    # plt.close()



