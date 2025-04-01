import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor




for re in ["promoters","enhancers"]:
    # read tfbs
    df_tfbs=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_{re}_k562.csv",index_col=0)
    # subset for proper background TFs
    proteins_background=pd.read_csv('MED_experiment_K562proteomics.txt', sep='\t',header=None)[0].tolist()
    df_tfbs=df_tfbs[df_tfbs['protein'].isin(proteins_background)].reset_index(drop=True)
    # calculate GC content
    seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")
    df_tfbs["motif"]=df_tfbs.apply(lambda row: seq_extractor.get_seq((row["chromosome"],row["start"],row["end"])), axis=1)
    df_tfbs["gc_content"]=df_tfbs.motif.apply(lambda x: (x.count("G")+x.count("C"))/len(x))
    # read mediator interactors
    df_med=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/PPI_MED/2024-08-07_MED-TF_interactions.txt",sep="\t")
    df_med=df_med[df_med["significant"]].reset_index(drop=True)
    # for each bait, add a column to df_tfbs, indicate whether the protein is a mediator interactor
    for bait in df_med["bait"].unique():
        logger.info(bait)
        tf_interactors=df_med[df_med["bait"]==bait]["gene"].tolist()
        df_tfbs[bait]=df_tfbs["protein"].apply(lambda x: x in tf_interactors)
    # add a column "Mediator": the OR of all bait columns
    bait_cols=df_med["bait"].unique()
    df_tfbs["Mediator"]=df_tfbs[bait_cols].apply(lambda x: x.any(), axis=1)
    # negate the Mediator column to get Non_interactors
    df_tfbs["Non_interactors"]=~df_tfbs["Mediator"]
    # melt
    value_vars=bait_cols.tolist()+["Mediator","Non_interactors"]
    df_melted = df_tfbs.melt(id_vars=["gc_content"], value_vars=value_vars, var_name="Bait", value_name="Is_Interactor")
    df_melted = df_melted[df_melted["Is_Interactor"]].reset_index(drop=True)
    # plot
    plt.figure(figsize=(3,2.5))
    sns.violinplot(x="Bait", y="gc_content", data=df_melted, cut=0, scale="width", linewidth=0.5, inner="quartile")
    plt.xticks(rotation=90,fontsize=5)
    plt.yticks(fontsize=5)
    plt.xlabel("Mediator Interactor",fontsize=7)
    plt.ylabel("GC content",fontsize=7)
    plt.title(f"GC content of motifs in {re}",fontsize=7)
    plt.tight_layout()
    plt.savefig(f"gc_content_{re}.pdf")
    plt.close()

