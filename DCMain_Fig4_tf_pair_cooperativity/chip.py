import pandas as pd
import numpy as np


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics
from tf_cooperativity import assign_cooperativity


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}_pe.csv")
    df=assign_cooperativity(df,0.3,0.7)
    df_chip_coloc=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_{cell_line}_summary.csv")
    df_chip_coloc.rename(columns={'tf1':'protein1','tf2':'protein2'},inplace=True)

    df=pd.merge(df,df_chip_coloc,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')
    df.rename(columns={'median_abs':'median_chip_peak_distance',
                       'tf_pair_count':'chip_overlap_count',
                       },inplace=True)
    df["log2_odds_ratio"]=np.log2(df["odds_ratio"])

    plot_violin_with_statistics(
        df=df,
        x_col="cooperativity",
        y_col="median_chip_peak_distance",
        x_label="Cooperativity",
        y_label="Median distance\nbetween ChIP peaks (bp)",
        title="Median distance vs Cooperativity index",
        output_file=f"chip_distance_vs_ci_{cell_line}.pdf"
    )

    plot_violin_with_statistics(
        df=df,
        x_col="cooperativity",
        y_col="log2_odds_ratio",
        x_label="Cooperativity",
        y_label="Log2 (colocalization odds ratio)",
        title="Colocalization tendency vs Cooperativity index",
        output_file=f"chip_colocalization_vs_ci_{cell_line}.pdf"
    )
    
    # conclusion: many nonfunctional bindings resulting in noninformative colocalization