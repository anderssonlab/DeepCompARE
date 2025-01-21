import pandas as pd
import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from seq_ops import SeqExtractor
from plotting import plot_violin_with_statistics


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



prefix='/isdata/alab/people/pcr980/DeepCompare/Pd5_motif_info/motif_info_thresh_500_'

for cell_line in ["hepg2","k562"]:
    
    df_promoters=pd.read_csv(f'{prefix}promoters_{cell_line}.csv')
    df_enhancers=pd.read_csv(f'{prefix}enhancers_{cell_line}.csv')

    df=pd.concat([df_promoters.assign(dataset="promoters"),
                df_enhancers.assign(dataset="enhancers")])


    seq_extractor=SeqExtractor("/isdata/alab/people/pcr980/Resource/hg38.fa")


    df["motif"]=df.apply(lambda row: seq_extractor.get_seq((row["chromosome"],row["start"],row["end"])), axis=1)
    df["gc_content"]=df.motif.apply(lambda x: (x.count("G")+x.count("C"))/len(x))


    # get redundant and codependent tfs
    tfs_redundant=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_redundant_{cell_line}.txt",header=None).iloc[:,0].tolist()
    tfs_codependent=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tfs_codependent_{cell_line}.txt",header=None).iloc[:,0].tolist()

    df["tf_type"]="other"
    df.loc[df["protein"].isin(tfs_redundant),"tf_type"]="redundant"
    df.loc[df["protein"].isin(tfs_codependent),"tf_type"]="codependent"
    df["tf_type"]=pd.Categorical(df["tf_type"],categories=["redundant","other","codependent"],ordered=True)



    plot_violin_with_statistics(
        df, 
        "tf_type",
        "gc_content",
        "TF Type",
        "GC percentage",
        None,
        f"gc_content_{cell_line}.pdf"
    )