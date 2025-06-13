import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns




import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from tf_cooperativity import assign_cooperativity



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


# 1. histogram of tf_pair ci
for cell_type in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_cooperativity_index_{cell_type}_pe.csv")
    df=assign_cooperativity(df,1,0.9,0.3,0.7)
    # # 1. histogram of tf_pair ci
    plt.figure(figsize=(2.3,2))
    # thin frame
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.histplot(df['cooperativity_index'], bins=40)
    # add vertical line at 0.3, 0.7,
    plt.axvline(x=0.3, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(x=0.7, color='gray', linestyle='--', linewidth=0.5)
    # add text at 0.15, 0.5, 0.85
    plt.text(0.02, 23, 'Redundant', fontsize=5, color='black')
    plt.text(0.35, 23, 'Intermediate', fontsize=5, color='black') # 1100 for tf pair, 20 for tf
    plt.text(0.75, 23, 'Synergistic', fontsize=5, color='black')
    plt.xlabel('TF pair synergy score',fontsize=7)
    plt.ylabel('Frequency',fontsize=7)
    plt.title(f'{cell_type}',fontsize=7)
    plt.xticks([0,0.3,0.7,1],fontsize=5)
    plt.yticks(fontsize=5)
    plt.ylim(0, 25)
    plt.tight_layout()
    plt.savefig(f'hist_tf_ci_{cell_type}_pe.pdf')
    plt.close()


    # # 2. heatmap of tf_pair ci
    # fig=plt.figure(figsize=(50,50))
    # df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
    # # remove rows with only NA
    # df = df.dropna(how='all')
    # # remove columns with only NA
    # df = df.dropna(axis=1, how='all')
    # gs = gridspec.GridSpec(1, 2, width_ratios=[0.05, 1])  # Adjust width_ratios to control size of the heatmap and color bar
    # # Create axes for the heatmap and color bar
    # cbar_ax = fig.add_subplot(gs[0, 0])  # Color bar on the left
    # heatmap_ax = fig.add_subplot(gs[0, 1])  # Heatmap on the right
    # sns.heatmap(df, cmap="coolwarm", vmin=0, vmax=1, ax=heatmap_ax, cbar_ax=cbar_ax)
    # plt.margin=0.01
    # plt.savefig(f'heatmap_tf_pair_ci_{cell_type}.pdf')
    # plt.close()



