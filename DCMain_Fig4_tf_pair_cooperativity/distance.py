import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,mannwhitneyu
from adjustText import adjust_text

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/Scripts_python")
from plotting import plot_violin_with_statistics
from tf_cooperativity import assign_cooperativity



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42



for cell_line in ["hepg2","k562"]:
    df=pd.read_csv(f"/isdata/alab/people/pcr980/DeepCompare/Pd7_TF_cooperativity/tf_pair_cooperativity_index_{cell_line}_pe.csv")
    df=assign_cooperativity(df,0.3,0.7)
    # set distance to plot
    df["distance_to_plot"]=df["nonlinear_distance"]
    # for linear pairs, use linear distance
    df.loc[df["cooperativity"]=="Linear","distance_to_plot"]=df.loc[df["cooperativity"]=="Linear","linear_distance"]
    # Analysis 1: distance vs cooperativity index
    plot_violin_with_statistics(
        df=df,
        x_col="cooperativity",
        y_col="distance_to_plot",
        x_label="Cooperativity",
        y_label="Median distance\nbetween TFBS pair (bp)",
        title=None,
        output_file=f"distance_vs_cooperativity_{cell_line}.pdf"
    )






    # Analysis 2: fix one TF, compare distance between codependent and redundant pairs

    # test_results = []
    # for protein in df["protein2"].unique():
    #     # Subset data for the current protein
    #     df_subset = df[df["protein2"] == protein].reset_index(drop=True)
    #     redundant_distances = df_subset.loc[df_subset["cooperativity"] == "redundant", "distance"]
    #     codependent_distances = df_subset.loc[df_subset["cooperativity"] == "codependent", "distance"]
    #     if len(redundant_distances) >= 2 and len(codependent_distances) >= 2:
    #         stat, p_value = mannwhitneyu(redundant_distances, codependent_distances, alternative="two-sided")
    #         test_results.append({"protein": protein, 
    #                             "statistic": stat, 
    #                             "p_value": p_value,
    #                             "redundant_median": np.median(redundant_distances),
    #                             "codependent_median": np.median(codependent_distances)})

    # # Convert test results into a DataFrame
    # df_res = pd.DataFrame(test_results)
    # df_res["significant"] = df_res["p_value"] < 0.05
    # df_res["significant"].sum()
    # # scatter plot
    # plt.figure(figsize=(6,6))
    # sns.scatterplot(data=df_res,x="redundant_median",y="codependent_median",hue="significant")
    # plt.xlabel("Median distance of redundant pairs")
    # plt.ylabel("Median distance of codependent pairs")
    # # add diagonal line
    # min_val=min(df_res[["codependent_median","redundant_median"]].min())
    # max_val=max(df_res[["codependent_median","redundant_median"]].max())
    # plt.plot([min_val,max_val],[min_val,max_val],color="black",linestyle="--")
    # # annotate each tf
    # texts = []
    # for i in range(df_res.shape[0]):
    #     if np.abs([df_res["redundant_median"][i]-df_res["codependent_median"][i]])>60:
    #         texts.append(plt.text(df_res["redundant_median"][i],df_res["codependent_median"][i],df_res["protein"][i]))

    # adjust_text(texts)
    # plt.savefig(f"distance_codependent_vs_redundant_{cell_line}.pdf")
    # plt.close()