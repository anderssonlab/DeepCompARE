import pandas as pd


import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics
from tf_cooperativity import assign_cooperativity




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


#---------------------------------
# Plot
#---------------------------------
cell_line = "hepg2"

for cell_line in ["hepg2", "k562"]:
    
    df_coop = pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}_pe.csv")
    df_coop =assign_cooperativity(df_coop,0.3,0.7)

    # protPair is combination of two protein names ordered alphabetically
    df_coop["protPair"]=df_coop.apply(lambda x: '_'.join(sorted([x["protein1"], x["protein2"]])), axis=1)
    # group by protPair and avg c_sum and cooperativity_index
    df_coop = df_coop.groupby("protPair").agg({"cooperativity_index": "mean", 
                                                "c_sum": "mean",
                                                "linearity_index": "mean"
                                                }).reset_index()
    
    # read ppi
    df_ppi = pd.read_csv(f"Pd1_Petra_data/2025-03-07_s10_PublishedPPIandProtComplexes_{cell_line}_pe.txt", sep='\t')
    df_ppi.rename(columns={'SWI/SNF': 'SWI_SNF'}, inplace=True)
    df = pd.merge(df_coop, df_ppi, on="protPair", how="inner")
    # Plot 1: Reported_PPI
    # make Reported_PPI categorical
    df["Reported_PPI"] = pd.Categorical(df["Reported_PPI"], categories=["No", "Yes"], ordered=True)
   
    plot_violin_with_statistics(
        # select non-nan cooperativity_index
        df=df[df["cooperativity_index"].notna()],
        x_col="Reported_PPI",
        y_col="cooperativity_index",
        x_label="Reported PPI",
        y_label="Cooperativity Index",
        title=None,
        output_file=f"ppi_reported_ppi_ci_{cell_line}.pdf"
    )
    
    plot_violin_with_statistics(
    df=df,
    x_col="Reported_PPI",
    y_col="linearity_index",
    x_label="Reported PPI",
    y_label="Linearity Index",
    title=None,
    output_file=f"ppi_reported_ppi_li_{cell_line}.pdf"
)

    # Plot 2: Protein complexes
    for col in ['Mediator', 'TFIIB', 'TFIID', 'TFIIE', 'TFIIF', 'TFIIH', 'SWI_SNF', 'POLII']:
        # make col categorical
        df[col] = pd.Categorical(df[col], categories=[0,1,2], ordered=True)
        plot_violin_with_statistics(
            df=df[df["cooperativity_index"].notna()],
            x_col=col,
            y_col="cooperativity_index",
            x_label=col,
            y_label="Cooperativity Index",
            title=None,
            output_file=f"ppi_{col}_ci_{cell_line}.pdf"
        )
        
        plot_violin_with_statistics(
            df=df,
            x_col=col,
            y_col="linearity_index",
            x_label=col,
            y_label="Linearity Index",
            title=None,
            output_file=f"ppi_{col}_li_{cell_line}.pdf"
        )
        
        

