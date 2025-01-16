import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#---------------------------------
# Read and merge data
#---------------------------------
cell_line = "hepg2"

for cell_line in ["hepg2", "k562"]:
    df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}.csv")
    # protPair is combination of two protein names ordered alphabetically
    df_coop["protPair"]=df_coop.apply(lambda x: '_'.join(sorted([x["protein1"], x["protein2"]])), axis=1)
    # group by protPair and avg c_sum and cooperativity_index
    df_coop = df_coop.groupby("protPair").agg({"cooperativity_index": "mean", "c_sum": "mean"}).reset_index()
    # select only c_sum > 1
    df_coop = df_coop[df_coop["c_sum"] > 1].reset_index(drop=True)
    # read ppi
    df_ppi = pd.read_csv(f"Pd1_Petra_data/2025-01-09_s10_PublishedPPIandProtComplexes_0.5_0.1_{cell_line}.txt", sep='\t')
    df_ppi.rename(columns={'SWI/SNF': 'SWI_SNF'}, inplace=True)
    df = pd.merge(df_coop, df_ppi, on="protPair", how="inner")
    # Plot 1: Reported_PPI
    plot_violin_with_statistics(
        df=df,
        x_col="Reported_PPI",
        y_col="cooperativity_index",
        x_label="Reported PPI",
        y_label="Cooperativity Index",
        title=None,
        output_file=f"ppi_reported_ppi_{cell_line}.pdf"
    )

    # Plot 2: Protein complexes
    for col in ['Mediator', 'TFIIB', 'TFIID', 'TFIIE', 'TFIIF', 'TFIIH', 'SWI_SNF', 'POLII']:
        plot_violin_with_statistics(
            df=df,
            x_col=col,
            y_col="cooperativity_index",
            x_label=col,
            y_label="Cooperativity Index",
            title=None,
            output_file=f"ppi_{col}_{cell_line}.pdf"
        )