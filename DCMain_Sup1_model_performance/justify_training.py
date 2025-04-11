import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt




import matplotlib
matplotlib.rcParams['pdf.fonttype']=42








df_mt=pd.read_csv("Pd2_benchmark_metrics/MT_metrics.csv")
# df_mt.file to categorical
df_mt['file'] = pd.Categorical(df_mt['file'], categories=sorted(df_mt['file'].unique()), ordered=True)



#------------------------------------------
# Helper functions
#------------------------------------------

def plot_performance(df, output_filename, title, color_col,
                     x_label="Dataset", y_label="Pearson correlation",
                     jitter_step=0.2):
    # Group by color_col and 'file' and compute mean, std, and count
    df_grouped = df.groupby([color_col, 'file'])['pcc'].agg(
        mean_pcc = 'mean',
        std_pcc  = 'std',
        count_pcc= 'count'
    ).reset_index()
    # Compute standard error and 95% confidence intervals
    df_grouped['sem_pcc'] = df_grouped['std_pcc'] / np.sqrt(df_grouped['count_pcc'])
    df_grouped['ci_low']  = df_grouped['mean_pcc'] - 1.96 * df_grouped['sem_pcc']
    df_grouped['ci_high'] = df_grouped['mean_pcc'] + 1.96 * df_grouped['sem_pcc']
    # Create unique x-positions for each file
    files = sorted(df_grouped['file'].unique())
    file_to_x = {f: i for i, f in enumerate(files)}
    # Get unique data categories and sort them for plotting order
    data_categories = sorted(df_grouped[color_col].unique())
    plt.figure(figsize=(3, 2.5))
    # Set thin frame and tick properties for better aesthetics
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis='both', which='both', width=0.5)
    # Initialize a jitter offset (will increment for each data category)
    jitter = 0
    colors=['tab:blue', 'tab:orange']
    for i, data_name in enumerate(data_categories):
        df_subset = df_grouped[df_grouped[color_col] == data_name]
        for _, row in df_subset.iterrows():
            x_val = file_to_x[row['file']] + jitter
            mean_y = row['mean_pcc']
            yerr_lower = mean_y - row['ci_low']
            yerr_upper = row['ci_high'] - mean_y
            plt.scatter(
                x_val,
                mean_y,
                color=colors[i],
                s=5,
                label=None  # prevent repeated labels
            )
            plt.errorbar(
                x_val,
                mean_y,
                yerr=[[yerr_lower], [yerr_upper]],
                fmt='none',
                ecolor=colors[i],
                capsize=0,
                elinewidth=1,
            )
        jitter += jitter_step
    # Create legend handles (one per category)
    for i,data_name in enumerate(data_categories):
        plt.scatter([], [], color=colors[i], label=data_name, s=10)
    #
    plt.legend(fontsize=5)
    # Configure x and y ticks
    plt.xticks(range(len(files)), files, fontsize=5, rotation=45)
    plt.yticks(fontsize=5)
    # Set axes labels and title
    plt.xlabel(x_label, fontsize=7, labelpad=0)
    plt.ylabel(y_label, fontsize=7)
    plt.title(title, fontsize=7)
    plt.tight_layout()
    # Save and close the figure
    plt.savefig(output_filename, dpi=300)
    plt.close()
    




#------------------------------------------
# Reg: clean data no worse than loose data
#------------------------------------------

# select model_type=="Reg"
df_reg=df_mt[df_mt['model_type']=='Reg'].copy().reset_index(drop=True)

# remove  df_reg.data_dir start with "Dataset_5cv_with_multilabel"
df_reg=df_reg[~df_reg['data_dir'].str.startswith("Dataset_5cv_with_multilabel")].copy().reset_index(drop=True)
df_reg['data']=df_reg['data_dir'].apply(lambda x: 'Credible regions' if x.startswith("Dataset_5cv_without_multilabel") else 'Loose positive')


plot_performance(df_reg, "credible_vs_loose.pdf", "Credible region","data")



#---------------------
# CR better than Reg
#----------------------
# select data_dir start with "Dataset_5cv_without_multilabel"
df_sub=df_mt[df_mt['data_dir'].str.startswith("Dataset_5cv_without_multilabel")].copy().reset_index(drop=True)
# remove model_type=="Class"
df_sub=df_sub[df_sub['model_type']!='Class'].copy().reset_index(drop=True)
# change model_type: CR-> "With classification", Reg->"Without classification"
df_sub['model_type'] = df_sub['model_type'].replace({'CR': 'With classification', 'Reg': 'Without classification'})

plot_performance(df_sub, "cr_vs_reg.pdf", "Classification as regulatization", "model_type")



#------------------------------------
# contrastive learning better than no contrastive learning
#------------------------------------

df_cr=df_mt[df_mt['model_type']=='CR'].copy().reset_index(drop=True)
# add column "data"
# if data_dir start with Dataset_5cv_without_multilabel, data='Without contrastive learning'
# if data_dir start with Dataset_5cv_with_multilabel, data='With contrastive learning'
df_cr['data'] = df_cr['data_dir'].apply(lambda x: 'With contrastive learning' if x.startswith("Dataset_5cv_with_multilabel") else 'Without contrastive learning')

plot_performance(df_cr, "contrastive_learning.pdf", "Contrast between cell types", "data")

