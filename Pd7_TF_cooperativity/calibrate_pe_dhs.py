import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, wilcoxon
import statsmodels.api as sm
import numpy as np



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42


# tf pair level
for cell_line in ['hepg2','k562']:
    df_pe=pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_pe.csv')
    # select c_sum>1
    df_pe=df_pe[df_pe['c_sum']>1].reset_index(drop=True)
    df_dhs=pd.read_csv(f'tf_pair_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs=df_dhs[df_dhs['c_sum']>1].reset_index(drop=True)
    df_merge=pd.merge(df_pe,df_dhs,on=['protein1','protein2'],
                    suffixes=('_pe','_dhs'),
                        how='inner')
    df_merge['c_sum_pe_clipped'] = df_merge['c_sum_pe'].clip(upper=50)
     # are the ranks preserved? Wilcoxon signed-rank test
    wilcoxon(df_merge['cooperativity_index_pe'],df_merge['cooperativity_index_dhs'])
    plt.figure(figsize=(2.3, 2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    hue="c_sum_pe_clipped",
                    s=5)
    plt.legend(fontsize=5,title='Sum(cooperativity)',title_fontsize=5,markerscale=0.5)
    plt.xlabel('CI defined by promoter + enhancer',fontsize=7)
    plt.ylabel('CI defined by DHS', fontsize=7)
    # annotate with pearson r
    r,p=pearsonr(df_merge['cooperativity_index_pe'],df_merge['cooperativity_index_dhs'])
    # add r at center of the plot
    plt.text(0.1,0.9,
             f'r={r:.2f}',
             horizontalalignment='center',
             verticalalignment='center',
             transform=plt.gca().transAxes,
             fontsize=5)
    # ticks font size 5
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_pair_{cell_line}.pdf')
    plt.close()





# Loop over cell lines
for cell_line in ['hepg2', 'k562']:
    # Read and filter the promoter+enhancer data
    df_pe = pd.read_csv(f'tf_cooperativity_index_{cell_line}_pe.csv')
    df_pe = df_pe[df_pe['c_sum'] > 1].reset_index(drop=True)
    # Read and filter the DHS data
    df_dhs = pd.read_csv(f'tf_cooperativity_index_{cell_line}_dhs.csv')
    df_dhs = df_dhs[df_dhs['c_sum'] > 1].reset_index(drop=True)
    # Merge datasets on protein2
    df_merge = pd.merge(df_pe, df_dhs, on=['protein2'],
                        suffixes=('_pe', '_dhs'), how='inner')
    wilcoxon(df_merge['cooperativity_index_pe'],df_merge['cooperativity_index_dhs'])
    # Clip c_sum_pe to an upper limit (here 50)
    df_merge['c_sum_pe'] = df_merge['c_sum_pe'].clip(upper=50)
    # Create the scatter plot with hue
    plt.figure(figsize=(2.3, 2.3))
    # thin frame
    plt.gca().spines['top'].set_linewidth(0.5)
    plt.gca().spines['right'].set_linewidth(0.5)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    sns.scatterplot(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    hue='c_sum_pe',
                    s=5)
    plt.xlabel('CI defined by promoter + enhancer',fontsize=7)
    plt.ylabel('CI defined by DHS', fontsize=7)
    plt.xticks(fontsize=5)
    plt.yticks(fontsize=5)
    # Fit a LOESS curve to the data. Adjust the "frac" parameter as needed.
    loess_result = sm.nonparametric.lowess(
        endog=df_merge['cooperativity_index_dhs'],
        exog=df_merge['cooperativity_index_pe'],
        frac=0.3
    )
    # Plot the LOESS curve
    plt.plot(loess_result[:, 0], loess_result[:, 1],
             color='tab:blue', linewidth=0.5)
    # For each desired x position, draw short dotted line segments (vertical and horizontal)
    # and note the LOESS predicted y values.
    x_positions = [0.3, 0.7]
    predicted_y = []  # We'll collect the predicted y values
    for x_val in x_positions:
        # Find the index in loess_result with x closest to x_val
        idx = np.abs(loess_result[:, 0] - x_val).argmin()
        y_val = loess_result[idx, 1]
        predicted_y.append(y_val)
        # Draw a vertical segment from (x_val, 0) to (x_val, y_val)
        plt.plot([x_val, x_val], [0, y_val], color='tab:gray', linestyle='--', linewidth=0.5)
        # Draw a horizontal segment from (0, y_val) to (x_val, y_val)
        plt.plot([0, x_val], [y_val, y_val], color='tab:gray', linestyle='--', linewidth=0.5)
    # Instead of annotating the crossing points with plt.text,
    # add the LOESS predicted y values to the y-axis ticks.
    ax = plt.gca()
    current_yticks = ax.get_yticks().tolist()
    # Append our predicted y values if they are not already present.
    for y_val in predicted_y:
        if y_val not in current_yticks:
            current_yticks.append(y_val)
    # Sort the ticks and set them back. Format tick labels with two decimals.
    current_yticks = np.sort(current_yticks)
    ax.set_yticks(current_yticks,fontsize=5)
    ax.set_yticklabels([f'{tick:.2f}' for tick in current_yticks],fontsize=5)
    ax.set_xticks([0,0.3,0.5,0.7,1])
    # Compute and annotate with Pearson r and p-value in the title
    r, p = pearsonr(df_merge['cooperativity_index_pe'],
                    df_merge['cooperativity_index_dhs'])
    plt.text(0.1, 0.9, f'r={r:.2f}',horizontalalignment='center',verticalalignment='center',transform=plt.gca().transAxes,fontsize=5)
    plt.xlim(0,df_merge['cooperativity_index_pe'].max()+0.1)
    plt.ylim(0,df_merge['cooperativity_index_dhs'].max()+0.1)
    plt.legend(fontsize=5,title='Sum(cooperativity)',title_fontsize=5, markerscale=0.5)
    # Add legend so the LOESS curve is identifiable
    # Save the figure
    plt.tight_layout()
    plt.savefig(f'calibrate_pe_dhs_{cell_line}.pdf')
    plt.close()



# Conclusion:DHS threshold:
# HepG2: 0.50,0.81
# K562: 0.46, 0.83



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
