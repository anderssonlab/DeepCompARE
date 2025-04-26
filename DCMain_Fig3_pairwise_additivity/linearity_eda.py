import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from loguru import logger

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#------------------------------------------------
# Pairwise independence score distribution and distance comparison
#------------------------------------------------


def plot_linear_nonlinear_distance(df,cell):
    df=df[df["cooperativity_index"].notnull()].reset_index(drop=True)
    df=df[df["linear_count"]>10].reset_index(drop=True)
    df=df[df["nonlinear_count"]>10].reset_index(drop=True)
    plt.figure(figsize=(2,2))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(data=df,x="linear_distance",y="nonlinear_distance",alpha=0.5, s=3)
    # count how many points are below the diagonal line
    logger.info(f"Percentage of points below diagonal line: {df[df['linear_distance']<df['nonlinear_distance']].shape[0]/df.shape[0]}")
    plt.xlabel("Independent cases", fontsize=7)
    plt.ylabel('Cooperative cases', fontsize=7)
    plt.title('Median distance (bp)', fontsize=7)
    # add diagonal line
    min_val=min(df["linear_distance"].min(),df["nonlinear_distance"].min())
    max_val=max(df["linear_distance"].max(),df["nonlinear_distance"].max())
    plt.plot([min_val,max_val],[min_val,max_val],color='grey',linestyle='--',linewidth=0.5)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f"distance_linear_vs_nonlinear_{cell}.pdf")
    plt.close() 



def plot_hist(df,cell,col,xlabel,outpath):
    plt.figure(figsize=(2.3,2))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.histplot(df[col],bins=50)
    # font size 5-7
    plt.xlabel(xlabel, fontsize=7)
    plt.ylabel('Frequency', fontsize=7)
    plt.title(cell, fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()




for cell in ['k562','hepg2']:
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell}_pe.csv')
    df["total_count"]=df["nonlinear_count"]+df["linear_count"]
    df=df[df["total_count"]>10].reset_index(drop=True)
    df_independence=df[df["linearity_index"]>0.9].reset_index(drop=True)
    logger.info(f"Independent TF pair {cell} percentage: {df_independence.shape[0]/df.shape[0]}")
    # plot distribution of TF pair linearity index
    plot_hist(df,cell,"linearity_index","Independence score",f'hist_tf_pair_li_{cell}.pdf')
    # for TFs with both linear and nonlinear counts, are linear distance larger than nonlinear distance?
    plot_linear_nonlinear_distance(df,cell)




