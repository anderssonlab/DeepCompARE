import pandas as pd

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from stat_tests import bin_and_label
from plotting import plot_violin_with_statistics


for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
    df=df[df['c_sum']>1].reset_index(drop=True)

    df_chip_coloc=pd.read_csv(f"/isdata/alab/people/pcr980/Resource/ReMap2022/Colocolization/colocolization_{cell_line}_summary.csv")
    df_chip_coloc.rename(columns={'tf1':'protein1','tf2':'protein2'},inplace=True)

    df=pd.merge(df,df_chip_coloc,left_on=['protein1','protein2'],right_on=['protein1','protein2'],how='inner')
    df.rename(columns={'mean_abs':'mean_chip_peak_distance',
                    'tf_pair_count':'chip_overlap_count',
                    },inplace=True)

    df=bin_and_label(df,"cooperativity_index",[0,0.3,0.7,1])
    plot_violin_with_statistics(
        df=df,
        x_col="cooperativity_index_bin",
        y_col="mean_chip_peak_distance",
        x_label="Cooperativity index bin",
        y_label="Mean distance between ChIP peaks (bp)",
        title="Mean distance vs Cooperativity Index",
        output_file=f"chip_distance_vs_ci_{cell_line}.pdf"
    )