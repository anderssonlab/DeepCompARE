
import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text



def plot_scatter_with_bias(data, x_column, y_column, x_label, y_label, title, output_filename,
                           x_gt_y_thresh, y_gt_x_thresh):
    # Compute bias
    data["bias_value"] = data[x_column] - data[y_column]

    # Plot setup
    plt.figure(figsize=(3.1, 2.8))
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)

    # Scatter plot with coolwarm colormap (no alpha)
    scatter = plt.scatter(
        x=data[x_column],
        y=data[y_column],
        c=data["bias_value"],
        cmap="coolwarm",
        s=10,
        edgecolor="none"
    )

    # Annotate significant points
    texts = []
    for _, row in data.iterrows():
        x = row[x_column]
        y = row[y_column]
        if x - y > x_gt_y_thresh or abs(x + y) > 0.95 or y - x > y_gt_x_thresh:
            texts.append(plt.text(x, y, row["protein"], fontsize=5, alpha=0.7))

    adjust_text(texts)

    # Axes and guidelines
    plt.xlabel(x_label, fontsize=7)
    plt.ylabel(y_label, fontsize=7)
    plt.title(title, fontsize=7)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.axvline(x=0, color="black", linewidth=0.5, linestyle="--")
    plt.axhline(y=0, color="black", linewidth=0.5, linestyle="--")

    # Colorbar setup
    cbar = plt.colorbar(scatter)
    cbar.set_label("Difference", fontsize=6)
    cbar.ax.tick_params(labelsize=6)

    # Save and close
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


# Load datasets
hepg2_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_hepg2_pe.csv'
k562_file = '/isdata/alab/people/pcr980/DeepCompare/Pd8_TF_effect_and_constraint/tf_effect_and_constraints_k562_pe.csv'

hepg2_data = pd.read_csv(hepg2_file)
k562_data = pd.read_csv(k562_file)

# Merge the files by "protein" using an inner join
merged_data = pd.merge(hepg2_data, k562_data, on="protein", suffixes=("_hepg2", "_k562"), how="inner")



plot_scatter_with_bias(
    data=merged_data,
    x_column="dstat_isa_cage_activity_hepg2",
    y_column="dstat_isa_cage_activity_k562",
    x_label="CAGE (HepG2)",
    y_label="CAGE (K562)",
    title="ISA-derived importance (D statistic)",
    output_filename="scatter_with_tf_names_hepg2_vs_k562.pdf",
    x_gt_y_thresh=0.45,
    y_gt_x_thresh=0.3
)



plot_scatter_with_bias(
    data=k562_data,
    x_column="dstat_isa_cage_activity",
    y_column="dstat_isa_starr_activity",
    x_label="CAGE (K562)",
    y_label="STARR (K562)",
    title="ISA-derived importance (D statistic)",
    output_filename="scatter_with_tf_names_cage_vs_starr.pdf",
    x_gt_y_thresh=0.3,
    y_gt_x_thresh=0.4
)




plot_scatter_with_bias(
    data=merged_data,
    x_column="dstat_isa_cage_activity_k562",
    y_column="dstat_isa_dhs_activity_k562",
    x_label="CAGE (K562)",
    y_label="DNase (K562)",
    title="ISA-derived importance (D statistic)",
    output_filename="scatter_with_tf_names_cage_vs_dhs.pdf",
    x_gt_y_thresh=0.25,
    y_gt_x_thresh=0.3
)



