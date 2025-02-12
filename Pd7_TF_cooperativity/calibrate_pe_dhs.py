import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import statsmodels.api as sm
import numpy as np



import matplotlib
matplotlib.rcParams['pdf.fonttype']=42

# TODO: how well are the ranks preserved?

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
    sns.scatterplot(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    hue="c_sum_pe_clipped",
                    s=5)
    plt.xlabel('CI defined by promoter + enhancer')
    plt.ylabel('CI defined by DHS')
    # annotate with pearson r
    r,p=pearsonr(df_merge['cooperativity_index_pe'],df_merge['cooperativity_index_dhs'])
    plt.title(f'TF pair level (r={r:.2f}, p={p:.2e})')
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
    # Clip c_sum_pe to an upper limit (here 50)
    df_merge['c_sum_pe'] = df_merge['c_sum_pe'].clip(upper=50)
    # Create the scatter plot with hue
    sns.scatterplot(data=df_merge,
                    x='cooperativity_index_pe',
                    y='cooperativity_index_dhs',
                    hue='c_sum_pe')
    plt.xlabel('CI defined by promoter + enhancer')
    plt.ylabel('CI defined by DHS')
    # Fit a LOESS curve to the data. Adjust the "frac" parameter as needed.
    loess_result = sm.nonparametric.lowess(
        endog=df_merge['cooperativity_index_dhs'],
        exog=df_merge['cooperativity_index_pe'],
        frac=0.3
    )
    # Plot the LOESS curve
    plt.plot(loess_result[:, 0], loess_result[:, 1],
             color='tab:blue', linewidth=2, label='LOESS')
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
        plt.plot([x_val, x_val], [0, y_val], color='tab:gray', linestyle='--')
        # Draw a horizontal segment from (0, y_val) to (x_val, y_val)
        plt.plot([0, x_val], [y_val, y_val], color='tab:gray', linestyle='--')
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
    ax.set_yticks(current_yticks)
    ax.set_yticklabels([f'{tick:.2f}' for tick in current_yticks])
    # Compute and annotate with Pearson r and p-value in the title
    r, p = pearsonr(df_merge['cooperativity_index_pe'],
                    df_merge['cooperativity_index_dhs'])
    plt.title(f'TF level (r={r:.2f}, p={p:.2e})')
    plt.xlim(0,df_merge['cooperativity_index_pe'].max()+0.1)
    plt.ylim(0,df_merge['cooperativity_index_dhs'].max()+0.1)
    # Add legend so the LOESS curve is identifiable
    plt.legend()
    # Save the figure
    plt.savefig(f'calibrate_pe_dhs_{cell_line}.pdf')
    plt.close()



# Conclusion:DHS threshold:
# HepG2: 0.50,0.81
# K562: 0.46, 0.83
