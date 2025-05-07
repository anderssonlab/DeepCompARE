import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import statsmodels.api as sm
import numpy as np


from matplotlib.colors import Normalize


import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


# tf pair level
for cell_line in ['hepg2','k562']:

    df_pe = pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_pe.csv')
    df_pe = df_pe[df_pe['c_sum'] > 1].reset_index(drop=True)
    df_dhs = pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs = df_dhs[df_dhs['c_sum'] > 1].reset_index(drop=True)
    df_merge = pd.merge(df_pe, df_dhs, on=['protein1', 'protein2'],
                        suffixes=('_pe', '_dhs'),
                        how='inner')
    df_merge['c_sum_pe_clipped'] = df_merge['c_sum_pe'].clip(upper=50)

    # Set up normalization and colormap
    norm = Normalize(vmin=10, vmax=50)
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    plt.figure(figsize=(2.5, 2.3))
    ax = plt.gca()
    # Thin frame
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    # Scatter plot with mapped colors
    sc = ax.scatter(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    c=df_merge['c_sum_pe_clipped'],
                    cmap=cmap,
                    norm=norm,
                    s=1,
                    alpha=0.3)

    # Colorbar
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label('Sum(cooperativity)', fontsize=5)
    cb.ax.tick_params(labelsize=5)

    # Labels and annotation
    plt.xlabel('Defined by promoter + enhancer', fontsize=7)
    plt.ylabel('Defined by DHS', fontsize=7)

    r, p = pearsonr(df_merge['cooperativity_index_pe'], df_merge['cooperativity_index_dhs'])
    plt.text(0.1, 0.9, f'r={r:.2f}', horizontalalignment='center',
            verticalalignment='center', transform=ax.transAxes, fontsize=5)
    plt.title("TF pair synergy score calibration", fontsize=7)
    plt.xticks([0,0.3,0.7,1],fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_pair_{cell_line}.pdf')
    plt.close()








# tf level
for cell_line in ['hepg2', 'k562']:
    # Read and filter the promoter+enhancer data
    df_pe = pd.read_csv(f'tf_cooperativity_index_{cell_line}_pe.csv')
    df_pe = df_pe[df_pe['c_sum'] > 5].reset_index(drop=True)
    # Read and filter the DHS data
    df_dhs = pd.read_csv(f'tf_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs = df_dhs[df_dhs['c_sum'] > 5].reset_index(drop=True)
    # Merge datasets on protein2
    df_merge = pd.merge(df_pe, df_dhs, on=['protein2'],
                        suffixes=('_pe', '_dhs'), how='inner')
    # Clip c_sum_pe to an upper limit (here 50)
    df_merge['c_sum_pe'] = df_merge['c_sum_pe'].clip(upper=50)
    # Normalize for color mapping
    norm = Normalize(vmin=10, vmax=50)
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    # Create the scatter plot with color bar
    plt.figure(figsize=(2.5, 2.3))
    ax = plt.gca()

    # Thin frame
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    sc = ax.scatter(df_merge['cooperativity_index_pe'],
                    df_merge['cooperativity_index_dhs'],
                    c=df_merge['c_sum_pe'],
                    cmap=cmap,
                    norm=norm,
                    s=1,
                    alpha=0.6)

    # Color bar instead of legend
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label('Sum(cooperativity)', fontsize=5)
    cb.ax.tick_params(labelsize=5)

    plt.xlabel('Defined by promoter + enhancer', fontsize=7)
    plt.ylabel('Defined by DHS', fontsize=7)
    plt.title("Synergy score calibration", fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)

    # LOESS smoothing
    loess_result = sm.nonparametric.lowess(
        endog=df_merge['cooperativity_index_dhs'],
        exog=df_merge['cooperativity_index_pe'],
        frac=0.3
    )
    plt.plot(loess_result[:, 0], loess_result[:, 1],
             color='tab:blue', linewidth=1)

    # Add guide lines and predicted y ticks
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

    # Pearson correlation
    r, p = pearsonr(df_merge['cooperativity_index_pe'], df_merge['cooperativity_index_dhs'])
    plt.text(0.1, 0.9, f'r={r:.2f}', horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=5)

    plt.xlim(0, df_merge['cooperativity_index_pe'].max() + 0.1)
    plt.ylim(0, df_merge['cooperativity_index_dhs'].max() + 0.1)

    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_{cell_line}.pdf')
    plt.close()








# Conclusion:DHS threshold:
# HepG2: 0.48,0.78
# K562: 0.43, 0.80



# jaccard similarity between redundant/codependent TFs


def jaccard_similarity_list(list_a, list_b):
    # Convert lists to sets
    set_a = set(list_a)
    set_b = set(list_b)
    # Calculate the intersection and union of the sets
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    print(intersection/len(set_a))
    # Return the Jaccard similarity, handle the case when union is empty
    return intersection / union if union != 0 else 0

for cell_line in ['hepg2','k562']:
    for coop in ['redundant','codependent']:
        tfs_pe=pd.read_csv(f'tfs_{coop}_{cell_line}_pe.txt',header=None)[0].tolist()
        tfs_dhs=pd.read_csv(f'tfs_{coop}_{cell_line}_dhs.txt',header=None)[0].tolist()
        print(f'{coop} {cell_line} Jaccard similarity: {jaccard_similarity_list(tfs_pe,tfs_dhs):.2f}')
