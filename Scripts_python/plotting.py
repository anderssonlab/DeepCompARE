import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu




def format_text(s):
    replace_map = {
        "distal_ti": "Distal, tissue-invariant",
        "proximal_ti": "Proximal, tissue-invariant",
        "distal_ts": "Distal, tissue-specific",
        "proximal_ts": "Proximal, tissue-specific",
        "distal_hepg2": "Distal (HepG2)",
        "distal_k562": "Distal (K562)",
        "proximal_hepg2": "Proximal (HepG2)",
        "proximal_k562": "Proximal (K562)",
        "hepg2": "HepG2",
        "k562": "K562",
        "cage": "CAGE",
        "dhs": "DNase",
        "DHS": "DNase",
        "starr": "STARR",
        "sure": "SuRE",
        "SURE": "SuRE",
        "enhancer": "Enhancer",
        "promoter": "Promoter",
    }
    # Apply exact matches first
    if s in replace_map:
        return replace_map[s]
    # Apply keyword-based replacements
    for old, new in replace_map.items():
        s = s.replace(old, new)
    # Fallback: general formatting
    return s.replace("_", " ")






def plot_violin_with_statistics(
    figsize,
    df, 
    x_col, 
    y_col, 
    x_label, 
    y_label, 
    title, 
    rotation,
    output_file
):
    """
    Plots a violin plot for a given dataset, annotates with Mann-Whitney U test p-values (for adjacent categories),
    and saves the plot to a file. The p-values are placed at a uniform height with flat, non-connecting lines,
    and slight jitter is applied to avoid overlapping annotations.
    
    Parameters:
        df (pd.DataFrame): Input dataframe.
        x_col (str): Column name for the x-axis categorical variable. (Must be a pandas Categorical)
        y_col (str): Column name for the y-axis continuous variable.
        x_label (str): Label for the x-axis.
        y_label (str): Label for the y-axis.
        title (str): Title for the plot.
        output_file (str): File path to save the plot.
    """
    # Apply text formatting
    x_label = format_text(x_label)
    y_label = format_text(y_label)
    title = format_text(title) if title is not None else None

    # Define base colors:
    white = "white"
    cool = "#1f77b4"    # blue
    gray = "#7f7f7f"    # gray
    warm = "#d62728"    # red

    # Extract ordered categories from the categorical column.
    bins = df[x_col].cat.categories.tolist()
    bin_counts = df[x_col].value_counts()

    # Choose palette automatically based on categories.
    if bins == ["Independent", "Redundant", "Intermediate", "Synergistic"]: 
        custom_palette = {
            bins[0]: white,
            bins[1]: cool,
            bins[2]: gray,
            bins[3]: warm
        }
    elif bins == ["No", "Yes"]:
        custom_palette = {
            "No": white,
            "Yes": white
        }
    elif bins == [0, 1, 2]:
        custom_palette = {
            0: white,
            1: white,
            2: white,
        }
    else:
        # Fallback: use Seaborn's default palette for the number of bins
        custom_palette = sns.color_palette("deep", len(bins))

    # Calculate p-values only for adjacent categories in the ordered list.
    p_values = []
    for i in range(len(bins) - 1):
        bin1, bin2 = bins[i], bins[i+1]
        group1 = df[df[x_col] == bin1][y_col]
        group2 = df[df[x_col] == bin2][y_col]
        stat, p_value = mannwhitneyu(group1, group2, alternative="two-sided")
        p_values.append({"bin1": bin1, "bin2": bin2, "p_value": p_value})
    p_values_df = pd.DataFrame(p_values)

    # Create the violin plot using the selected palette.
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

    sns.violinplot(
        data=df, 
        x=x_col, 
        y=y_col, 
        order=bins, 
        cut=0, 
        linewidth=0.5, 
        alpha=0.5,
        inner="quart",
        palette=custom_palette,
        linecolor='black'
    )

    # Determine y-axis range and set a base annotation height
    y_min, y_max = df[y_col].min(), df[y_col].max()
    y_margin = (y_max - y_min) * 0.05
    annotation_base = y_max + y_margin

    # Annotate plot with adjacent p-values.
    np.random.seed(42)
    jitters = np.random.uniform(-y_margin * 0.2, y_margin * 0.2, len(p_values_df))

    for i, row in p_values_df.iterrows():
        bin1, bin2 = row["bin1"], row["bin2"]
        p_value = row["p_value"]
        x1, x2 = bins.index(bin1), bins.index(bin2)
        line_y = annotation_base + jitters[i]
        plt.plot([x1, x2], [line_y, line_y], lw=0.5, color="black")
        plt.text((x1 + x2) / 2, line_y + y_margin * 0.15, 
                 f"p={p_value:.1e}", ha="center", fontsize=5)

    # Add counts to x-axis labels.
    bin_labels = [f"{bin_}\n(n={bin_counts[bin_]})" for bin_ in bins]
    plt.xticks(ticks=range(len(bins)), labels=bin_labels, fontsize=5, rotation=rotation)
    plt.yticks(fontsize=5)

    # Finalize the plot with labels and layout adjustments.
    plt.xlabel(x_label, fontsize=7)
    plt.ylabel(y_label, fontsize=7)
    if title is not None:
        plt.title(title, fontsize=5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
