# plotting_utils.py

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import sys

sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python")
from plotting import format_text


def annotate_significance(ax, df, x_col, neighbor_pairs, showfliers, offset_frac=0.25):
    for (cat_a, cat_b) in neighbor_pairs:
        group_a = df.loc[df["region_type"] == cat_a, x_col].dropna()
        group_b = df.loc[df["region_type"] == cat_b, x_col].dropna()
        if len(group_a) == 0 or len(group_b) == 0:
            continue

        stat, pval = mannwhitneyu(group_a, group_b, alternative="two-sided")
        if pval > 0.05:
            continue

        if pval <= 1e-3:
            stars = "***"
        elif pval <= 1e-2:
            stars = "**"
        else:
            stars = "*"

        if showfliers:
            x_ref  = max(group_a.max(), group_b.max())
            x_line = x_ref * 1.01
            x_text = x_ref * 1.02
        else:
            pct90_a = group_a.quantile(0.90)
            pct90_b = group_b.quantile(0.90)
            x_ref   = max(pct90_a, pct90_b)
            x_line  = x_ref * 1.02
            x_text  = x_ref * 1.03

        y1 = df["region_type"].cat.categories.get_loc(cat_a)
        y2 = df["region_type"].cat.categories.get_loc(cat_b)
        dy = y2 - y1
        start_y = y1 + dy * offset_frac
        end_y   = y2 - dy * offset_frac

        ax.plot([x_line, x_line], [start_y, end_y], lw=0.2, color="black")
        ax.text(
            x_text,
            (y1 + y2) / 2,
            stars,
            va="center",
            ha="left",
            fontsize=5,
            fontweight="bold",
            color="red"
        )


def plot_region_boxplot(
    df,
    x_col,
    title,
    outfile,
    neighbor_pairs,
    color=None,
    box_alpha=0.3,
    showfliers=True,
    ax=None,
    figsize=(3, 1.4),
    xlim=None
):
    df = df.copy()
    df["region_type"] = df["region_type"].cat.rename_categories(lambda s: format_text(s))
    formatted_pairs = [(format_text(a), format_text(b)) for (a, b) in neighbor_pairs]

    own_fig = False
    if ax is None:
        own_fig = True
        plt.figure(figsize=figsize)
        ax = plt.gca()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)

    if color:
        boxprops = {"facecolor": color, "alpha": box_alpha}
        medianprops = {"color": color}
    else:
        boxprops = {"facecolor": "none"}
        medianprops = {"color": "black"}

    sns.boxplot(
        x=x_col,
        y="region_type",
        data=df,
        ax=ax,
        linewidth=0.5,
        boxprops=boxprops,
        whiskerprops={"color": "black"},
        showcaps=False,
        medianprops=medianprops,
        showfliers=showfliers,
        flierprops={"marker": "o", "markersize": 0.8},
    )

    annotate_significance(ax, df, x_col, formatted_pairs, showfliers)

    ax.set_ylabel("")
    if "count" in x_col:
        ax.set_xlabel("TFBS count", fontsize=6)
    else:
        ax.set_xlabel("Synergy score", fontsize=6)
    ax.set_title(title, fontsize=6)
    ax.tick_params(axis="x", labelsize=5)
    ax.tick_params(axis="y", labelsize=5)

    if xlim:
        ax.set_xlim(xlim)

    if own_fig:
        plt.tight_layout()
        plt.savefig(outfile)
        plt.close()
        print(f"Saved → {outfile}")


def plot_combined(
    panel_specs,
    neighbor_pairs,
    outfile,
    figsize=(5, 2.8),
    nrows=1,
    ncols=None,
    sharey=True,
    wspace=0.5,
    hspace=0.5
):
    total = len(panel_specs)
    if ncols is None:
        ncols = int((total + nrows - 1) // nrows)

    sns.set_style("white", {"axes.grid": False})
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize,
        sharey=sharey
    )
    plt.subplots_adjust(wspace=wspace, hspace=hspace)

    axes_list = [axes] if nrows * ncols == 1 else axes.flatten()

    for ax, spec in zip(axes_list, panel_specs):
        plot_region_boxplot(
            df             = spec["df"],
            x_col          = spec["x_col"],
            title          = spec["title"],
            outfile        = "unused.pdf",
            neighbor_pairs = neighbor_pairs,
            color          = spec.get("color", None),
            box_alpha      = spec.get("box_alpha", 0.3),
            showfliers     = spec.get("showfliers", True),
            ax             = ax,
            xlim           = spec.get("xlim", None)
        )

    for ax in axes_list[len(panel_specs):]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close(fig)
    print(f"Saved → {outfile}")
