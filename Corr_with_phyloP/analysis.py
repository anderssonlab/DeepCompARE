import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import read_featimp,remove_nan_inf
from stat_tests import bin_and_label, calc_odds_ratio



# use ism by default

#----------------------
# Plot odds ratio of significantly conserved bases
#----------------------


def read_ism_and_cons(file_prefix,phylop_prefix,track_num):
    # read in both files and merge into one dataframe
    df_ism=read_featimp(f"/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/ism_{file_prefix}.csv",track_num)
    df_cons=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Corr_with_phyloP/{phylop_prefix}_{file_prefix}.csv",header=None)
    ism=df_ism.values.reshape(-1,1).squeeze()
    cons=df_cons.values.reshape(-1,1).squeeze()
    assert ism.shape==cons.shape
    ism,cons=remove_nan_inf(ism,cons)
    df = pd.DataFrame({'ism':ism,'cons':cons})
    # split df to conserved and accelerated
    df_conserved=df[df["cons"]>0].copy().reset_index(drop=True)
    df_conserved=bin_and_label(df_conserved, 'cons', [0,2.27,np.inf])
    df_conserved["Evolution_type"]=df_conserved.replace({"Bin":{"0-2.27":"Non-sig","2.27-inf":"Sig"}})["Bin"]
    df_accelerated=df[df["cons"]<0].copy().reset_index(drop=True)
    df_accelerated=bin_and_label(df_accelerated, 'cons', [-np.inf,-2.27,0])
    df_accelerated["Evolution_type"]=df_accelerated.replace({"Bin":{"-inf--2.27":"Sig","-2.27-0":"Non-sig"}})["Bin"]
    return df_conserved,df_accelerated

def bin_ism_calc_OR(df,phylop_label,status_label):
    df.drop(columns=["Bin"],inplace=True) 
    df=bin_and_label(df, 'ism', [-np.inf,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,np.inf])
    df_plot=df.groupby(["Bin","Evolution_type"]).size().unstack()
    df_plot=df_plot[df_plot.sum(axis=1)>10]
    df_plot["or_Sig"]=df_plot.apply(lambda x: calc_odds_ratio(df_plot,x.name,"Sig"),axis=1)
    df_plot["phyloP"]=phylop_label
    df_plot["status"]=status_label
    df_plot.drop(columns=["Non-sig","Sig"],inplace=True)
    return df_plot


def analyze_re(file_prefix,track_num):

    df_100way_cons, df_100way_acc=read_ism_and_cons(file_prefix,"100way",track_num)
    df_241way_cons, df_241way_acc=read_ism_and_cons(file_prefix,"241way",track_num)
    df_447way_cons, df_447way_acc=read_ism_and_cons(file_prefix,"447way",track_num)
    _, df_447wayLRT_acc=read_ism_and_cons(file_prefix,"447wayLRT",track_num)

    or_100way_cons=bin_ism_calc_OR(df_100way_cons,"100way","conserved")
    or_241way_cons=bin_ism_calc_OR(df_241way_cons,"241way","conserved")
    or_447way_cons=bin_ism_calc_OR(df_447way_cons,"447way","conserved")

    or_100way_acc=bin_ism_calc_OR(df_100way_acc,"100way","accelerated")
    or_241way_acc=bin_ism_calc_OR(df_241way_acc,"241way","accelerated")
    or_447way_acc=bin_ism_calc_OR(df_447way_acc,"447way","accelerated")
    or_447wayLRT_acc=bin_ism_calc_OR(df_447wayLRT_acc,"447wayLRT","accelerated")

    # concatenate the above by row
    df_or=pd.concat([or_100way_cons,or_241way_cons,or_447way_cons,
                    or_100way_acc,or_241way_acc,or_447way_acc,or_447wayLRT_acc],
                    axis=0)
    df_or["Predicted effect size"]=df_or.index
    # negate or_Sig if status is accelerated
    df_or["or_Sig"]=df_or.apply(lambda x: -x["or_Sig"] if x["status"]=="accelerated" else x["or_Sig"],axis=1)

    plt.figure(figsize=(6, 8))
    plt.subplots_adjust(left=0.2, bottom=0.2)  # Increase left and lower margins
    sns.scatterplot(data=df_or[df_or["status"]=="conserved"],x="Predicted effect size", y="or_Sig", hue="phyloP",legend=False)
    sns.lineplot(data=df_or[df_or["status"]=="conserved"],x="Predicted effect size", y="or_Sig", hue="phyloP", linestyle="--",legend=False)
    sns.scatterplot(data=df_or[df_or["status"]=="accelerated"],x="Predicted effect size", y="or_Sig", hue="phyloP")
    sns.lineplot(data=df_or[df_or["status"]=="accelerated"],x="Predicted effect size", y="or_Sig", hue="phyloP", linestyle="--",legend=False)
    plt.xticks(rotation=45)
    # put legend in the middle
    plt.legend(loc='center right')  
    # add horizontal line: y=0
    plt.axhline(0, color='black', linestyle='--')
    plt.ylabel("Odds ratio of significant base")
    plt.title(f"Odds ratio in {file_prefix}, track {track_num}")
    plt.savefig(f"OR_{file_prefix}_{track_num}.png")
    plt.close()


if __name__=="__main__":
    logger.add("analysis.log")
    for track_num in [8,9,10,11,12,13,14,15]:
        for file_prefix in ["promoters_hepg2","enhancers_hepg2","promoters_k562","enhancers_k562"]:
            analyze_re(file_prefix,track_num)
            logger.info(f"Done with {file_prefix}, track {track_num}")
    



# nohup python3 analysis.py &

    


