import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec



import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity


# 1. histogram of tf_pair ci
for cell_type in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_type}_pe.csv")
    df=assign_cooperativity(df,1,0.9,0.3,0.7)
    # # 1. histogram of tf_pair ci
    # plt.figure(figsize=(2.3,2))
    # # thin frame
    # plt.gca().spines['top'].set_linewidth(0.5)
    # plt.gca().spines['bottom'].set_linewidth(0.5)
    # plt.gca().spines['left'].set_linewidth(0.5)
    # plt.gca().spines['right'].set_linewidth(0.5)
    # sns.histplot(df['cooperativity_index'], bins=40)
    # plt.xlabel('TF pair synergy score',fontsize=7)
    # plt.ylabel('Frequency',fontsize=7)
    # plt.xticks(fontsize=5)
    # plt.yticks(fontsize=5)
    # plt.tight_layout()
    # plt.savefig(f'plot_hist_tf_pair_ci_{cell_type}_pe.pdf')
    # plt.close()


    # 2. heatmap of tf_pair ci
    fig=plt.figure(figsize=(50,50))
    df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
    # remove rows with only NA
    df = df.dropna(how='all')
    # remove columns with only NA
    df = df.dropna(axis=1, how='all')
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.05, 1])  # Adjust width_ratios to control size of the heatmap and color bar
    # Create axes for the heatmap and color bar
    cbar_ax = fig.add_subplot(gs[0, 0])  # Color bar on the left
    heatmap_ax = fig.add_subplot(gs[0, 1])  # Heatmap on the right
    sns.heatmap(df, cmap="coolwarm", vmin=0, vmax=1, ax=heatmap_ax, cbar_ax=cbar_ax)
    plt.margin=0.01
    plt.savefig(f'plot_heatmap_tf_pair_ci_{cell_type}.pdf')
    plt.close()



