import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('/isdata/alab/people/pcr980/Scripts_python')
from tf_cooperativity import assign_cooperativity



def plot_hist(df,cell,col,xlabel,outpath):
    plt.figure(figsize=(3,2.5))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.histplot(df[col],bins=50)
    # font size 5-7
    plt.xlabel(xlabel, fontsize=7)
    plt.ylabel('Frequency', fontsize=7)
    plt.title(f'Linearity Index Distribution ({cell})', fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()
    


def plot_linear_nonlinear_distance(df,cell):
    df=df[df["cooperativity_index"].notnull()].reset_index(drop=True)
    df=df[df["linear_count"]>10].reset_index(drop=True)
    df=df[df["nonlinear_count"]>10].reset_index(drop=True)
    plt.figure(figsize=(2.5,2.5))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(data=df,x="linear_distance",y="nonlinear_distance",alpha=0.5, s=3)
    plt.xlabel("Median distance of linear cases", fontsize=7)
    plt.ylabel('Median distance of nonlinear cases', fontsize=7)
    plt.title('Distance between TFBS pairs', fontsize=7)
    # add diagonal line
    min_val=min(df["linear_distance"].min(),df["nonlinear_distance"].min())
    max_val=max(df["linear_distance"].max(),df["nonlinear_distance"].max())
    plt.plot([min_val,max_val],[min_val,max_val],color='grey',linestyle='--',linewidth=0.5)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f"distance_linear_vs_nonlinear_{cell}.pdf")



#-------------------------
# TF pair linearity
#-------------------------
for cell in ['k562','hepg2']:
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell}_pe.csv')
    df["total_count"]=df["nonlinear_count"]+df["linear_count"]
    df=df[df["total_count"]>10].reset_index(drop=True)
    # plot distribution of TF pair linearity index
    plot_hist(df,cell,"linearity_index","Linearity index",f'hist_tf_pair_li_{cell}.pdf')
    # for TFs with both linear and nonlinear counts, are linear distance larger than nonlinear distance?
    plot_linear_nonlinear_distance(df,cell)







#-------------------------
# Do linearity differ between enhancer and promoter?
#--------------------------

for cell in ['hepg2','k562']:
    df_promoter=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell}_promoter_pe.csv')
    df_promoter=assign_cooperativity(df_promoter,0.3,0.7)
    df_enhancer=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell}_enhancer_pe.csv')
    df_enhancer=assign_cooperativity(df_enhancer,0.3,0.7)
    # inner merge by TF pair
    df=pd.merge(df_promoter,df_enhancer,on=["protein1","protein2"],suffixes=('_promoter','_enhancer'))
    # scatter linearity_index_enhancer and linearity_index_promoter
    plt.figure(figsize=(2.5,2.5))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(data=df,x="linearity_index_promoter",y="linearity_index_enhancer",alpha=0.5, s=3)
    plt.xlabel("Linearity index (promoter)", fontsize=7)
    plt.ylabel('Linearity index (enhancer)', fontsize=7)
    plt.title('Linearity index comparison', fontsize=7)
    # add diagonal line
    min_val=min(df["linearity_index_promoter"].min(),df["linearity_index_enhancer"].min())
    max_val=max(df["linearity_index_promoter"].max(),df["linearity_index_enhancer"].max())
    plt.plot([min_val,max_val],[min_val,max_val],color='grey',linestyle='--',linewidth=0.5)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f"linearity_index_promoter_vs_enhancer_{cell}.pdf")
    plt.close()   


