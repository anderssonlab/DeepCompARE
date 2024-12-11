import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

for cell_type in ["hepg2","k562"]:
    # 1. histogram of tf_pair ci
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_type}.csv")
    df=df[df["c_sum"]>1].reset_index(drop=True)
    sns.histplot(df['cooperativity_index'], bins=100)
    plt.xlabel('TF pair cooperativity index')
    plt.ylabel('Frequency')
    plt.title('TF pair cooperativity index distribution')
    plt.savefig(f'plot_hist_tf_pair_ci_{cell_type}.pdf')
    plt.close()


    # 2. heatmap of tf_pair ci
    plt.figure(figsize=(50,50))
    df=df.pivot(index="protein1",columns="protein2",values="cooperativity_index")
    sns.heatmap(df,cmap="coolwarm",vmin=0,vmax=1)
    plt.title("TF pair cooperativity ratio")
    plt.subplots_adjust(bottom=0.3)
    plt.subplots_adjust(left=0.3)
    plt.savefig(f"plot_heatmap_tf_pair_ci_{cell_type}.pdf")
    plt.close()


