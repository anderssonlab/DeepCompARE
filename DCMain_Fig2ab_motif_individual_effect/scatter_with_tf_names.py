import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42




def plot_scatter_with_bias(data, x_column, y_column, output_filename):
    """
    Plots a scatter plot showing bias between two datasets based on specified columns.
    Parameters:
        data (pd.DataFrame): The dataset containing the columns for comparison.
        x_column (str): The column name for the x-axis values.
        y_column (str): The column name for the y-axis values.
        output_filename (str): The filename to save the resulting plot.
    """
    # Create generalized bias and transparency columns
    data["bias"] = data.apply(
        lambda row: f"{x_column}-biased" if row[x_column] > row[y_column]
                    else (f"{y_column}-biased" if row[x_column] < row[y_column]
                          else "Neutral"),
        axis=1
    )
    #
    data["alpha"] = data.apply(
        lambda row: min(1.0, max(0.1, 2*abs(row[x_column] - row[y_column]))),
        axis=1
    )
    #
    # Plot setup
    plt.figure(figsize=(2.8, 2.8))
    ax = plt.gca()
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    #
    scatter = plt.scatter(
        x=data[x_column],
        y=data[y_column],
        c=data["bias"].map({f"{x_column}-biased": "#D98880", f"{y_column}-biased": "#1E90FF", "Neutral": "grey"}), 
        alpha=data["alpha"],
        s=10,
        edgecolor="none"
    )
    #
    # Annotate significant points
    texts = []
    for _, row in data.iterrows():
        x = row[x_column]
        y = row[y_column]
        if abs(x - y) > 0.4 or abs(x + y) > 0.95:  # Opposite signs or large sum of coordinates
            texts.append(plt.text(x, y, row["protein"], fontsize=5, alpha=0.7))
    #
    # Adjust text (if `adjust_text` is available)
    adjust_text(texts)
    #
    # Axes labels and guidelines
    plt.xlabel(f"Effect on {x_column}", fontsize=7)
    plt.ylabel(f"Effect on {y_column}", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.axvline(x=0, color="black", linewidth=0.5, linestyle="--")
    plt.axhline(y=0, color="black", linewidth=0.5, linestyle="--")
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()




# Load datasets
hepg2_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_hepg2.csv'
k562_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_k562.csv'

hepg2_data = pd.read_csv(hepg2_file)
k562_data = pd.read_csv(k562_file)

# Merge the files by "protein" using an inner join
merged_data = pd.merge(hepg2_data, k562_data, on="protein", suffixes=("_hepg2", "_k562"), how="inner")


plot_scatter_with_bias(
    data=merged_data,
    x_column="dstat_isa_cage_activity_hepg2",
    y_column="dstat_isa_cage_activity_k562",
    output_filename="scatter_with_tf_names_hepg2_vs_k562.pdf"
)




plot_scatter_with_bias(
    data=k562_data,
    x_column="dstat_isa_cage_activity",
    y_column="dstat_isa_starr_activity",
    output_filename="scatter_with_tf_names_cage_vs_starr.pdf"
)

