import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import statsmodels.api as sm
import numpy as np
from matplotlib.colors import Normalize
import matplotlib
import sys

sys.path.insert(1, "/isdata/alab/people/pcr980/Scripts_python")
from plotting import format_text  # This is assumed to contain your format_text() function

matplotlib.rcParams['pdf.fonttype'] = 42

# Shared normalization and colormap for consistent colorbars
norm = Normalize(vmin=10, vmax=50)
cmap = sns.color_palette("rocket_r", as_cmap=True)

def preprocess(df, c_sum_thresh):
    df["total_count"] = df["linear_count"] + df["nonlinear_count"]
    df = df[df["total_count"] > 10].reset_index(drop=True)
    df = df[df['i_sum'] > c_sum_thresh].reset_index(drop=True)
    return df

# TF pair level
for cell_line in ['hepg2', 'k562']:
    label = format_text(cell_line)
    
    df_pe = pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_pe.csv')
    df_pe = preprocess(df_pe, c_sum_thresh=1)
    df_dhs = pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs = preprocess(df_dhs, c_sum_thresh=1)

    df_merge = pd.merge(df_pe, df_dhs, on=['protein1', 'protein2'],
                        suffixes=('_pe', '_dhs'), how='inner')
    df_merge['c_sum_pe_clipped'] = df_merge['c_sum_pe'].clip(upper=50)

    plt.figure(figsize=(2.6, 2.3))
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    sc = ax.scatter(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    c=df_merge['c_sum_pe_clipped'],
                    cmap=cmap,
                    norm=norm,
                    s=1,
                    alpha=0.3)

    cb = plt.colorbar(sc, ax=ax, ticks=[10, 20, 30, 40, 50])
    cb.set_label('Sum(cooperativity) (clipped at 50)', fontsize=5)
    cb.ax.tick_params(labelsize=5)

    plt.xlabel('Defined by promoter + enhancer', fontsize=7)
    plt.ylabel('Defined by DHS', fontsize=7)

    r, _ = pearsonr(df_merge['cooperativity_index_pe'], df_merge['cooperativity_index_dhs'])
    plt.text(0.1, 0.9, f'r={r:.2f}', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=5)

    plt.title(f"TF pair synergy score ({label})", fontsize=7)
    plt.xticks([0, 0.3, 0.7, 1], fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_pair_{cell_line}.pdf')
    plt.close()




# TF level
for cell_line in ['hepg2', 'k562']:
    label = format_text(cell_line)
    
    df_pe = pd.read_csv(f'tf_cooperativity_index_{cell_line}_pe.csv')
    df_pe = preprocess(df_pe, c_sum_thresh=5)
    df_dhs = pd.read_csv(f'tf_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs = preprocess(df_dhs, c_sum_thresh=5)

    df_merge = pd.merge(df_pe, df_dhs, on=['protein2'],
                        suffixes=('_pe', '_dhs'), how='inner')
    df_merge['c_sum_pe'] = df_merge['c_sum_pe'].clip(upper=50)

    plt.figure(figsize=(2.6, 2.3)) # (2.2, 2) for main
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    sc = ax.scatter(df_merge['cooperativity_index_pe'],
                    df_merge['cooperativity_index_dhs'],
                    c=df_merge['c_sum_pe'],
                    cmap=cmap,
                    norm=norm,
                    s=1,
                    alpha=0.6)

    cb = plt.colorbar(sc, ax=ax, ticks=[10, 20, 30, 40, 50])
    cb.set_label('Sum(cooperativity) (clipped at 50)', fontsize=5)
    cb.ax.tick_params(labelsize=5)

    plt.xlabel('Defined by promoter + enhancer', fontsize=7)
    plt.ylabel('Defined by DHS', fontsize=7)
    plt.title(f"TF synergy score ({label})", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)

    loess_result = sm.nonparametric.lowess(
        endog=df_merge['cooperativity_index_dhs'],
        exog=df_merge['cooperativity_index_pe'],
        frac=0.3
    )
    plt.plot(loess_result[:, 0], loess_result[:, 1], color='tab:blue', linewidth=1)

    x_positions = [0.3, 0.7]
    predicted_y = []
    for x_val in x_positions:
        idx = np.abs(loess_result[:, 0] - x_val).argmin()
        y_val = loess_result[idx, 1]
        predicted_y.append(y_val)
        plt.plot([x_val, x_val], [0, y_val], color='tab:gray', linestyle='--', linewidth=0.5)
        plt.plot([0, x_val], [y_val, y_val], color='tab:gray', linestyle='--', linewidth=0.5)

    current_yticks = ax.get_yticks().tolist()
    for y_val in predicted_y:
        if y_val not in current_yticks:
            current_yticks.append(y_val)
    current_yticks = np.sort(current_yticks)
    ax.set_yticks(current_yticks)
    ax.set_yticklabels([f'{tick:.2f}' for tick in current_yticks], fontsize=5)
    ax.set_xticks([0, 0.3, 0.7, 1])

    r, _ = pearsonr(df_merge['cooperativity_index_pe'], df_merge['cooperativity_index_dhs'])
    plt.text(0.1, 0.9, f'r={r:.2f}', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=5)

    plt.xlim(0, df_merge['cooperativity_index_pe'].max() + 0.1)
    plt.ylim(0, df_merge['cooperativity_index_dhs'].max() + 0.1)
    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_{cell_line}.pdf')
    plt.close()








# Jaccard similarity between redundant/codependent TFs
def jaccard_similarity_list(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    print(intersection / len(set_a))
    return intersection / union if union != 0 else 0

for cell_line in ['hepg2', 'k562']:
    for coop in ['redundant', 'codependent']:
        tfs_pe = pd.read_csv(f'tfs_{coop}_{cell_line}_pe.txt', header=None)[0].tolist()
        tfs_dhs = pd.read_csv(f'tfs_{coop}_{cell_line}_dhs.txt', header=None)[0].tolist()
        print(f'{coop} {cell_line} Jaccard similarity: {jaccard_similarity_list(tfs_pe, tfs_dhs):.2f}')
