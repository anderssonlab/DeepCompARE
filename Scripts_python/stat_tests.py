import pandas as pd
import numpy as np
from loguru import logger
from scipy.stats import mannwhitneyu
from utils import remove_nan_inf
from scipy.stats import pearsonr, fisher_exact
import matplotlib.pyplot as plt



def match_by_decile(df,score_col,label_col):
    logger.info(f"Original data frame has shape {df.shape}")
    df=df.copy()    
    # test if the score distribution is different between groups
    scores_true = df[df[label_col] == True][score_col]
    scores_false = df[df[label_col] == False][score_col]
    if len(scores_true) == 0 or len(scores_false) == 0:
        logger.warning("No positive or negative samples in the data frame")
        return pd.DataFrame()
    _, p_value = mannwhitneyu(scores_true, scores_false)
    if p_value > 0.1:
        logger.info(f"No significant difference in score distribution between groups, p-value={p_value:.3f}")
        return df
    df['score_decile'] = pd.qcut(df[score_col], 20, labels=False, duplicates='drop')
    
    balanced_dfs = []
    for decile in range(10):
        decile_df = df[df['score_decile'] == decile]
        min_size = decile_df[label_col].value_counts().min()
        balanced_df_decile = decile_df.groupby(label_col).apply(lambda x: x.sample(n=min_size)).reset_index(drop=True)
        balanced_dfs.append(balanced_df_decile)
    
    balanced_df = pd.concat(balanced_dfs).reset_index(drop=True)
    logger.info(f"Balanced data frame has shape {balanced_df.shape}")
    
    scores_true = balanced_df[balanced_df[label_col] == True][score_col]
    scores_false = balanced_df[balanced_df[label_col] == False][score_col]
    if len(scores_true) == 0 or len(scores_false) == 0:
        logger.warning("No positive or negative samples in the balanced data frame")
        return pd.DataFrame()
    
    _, p_value = mannwhitneyu(scores_true, scores_false)
    if p_value > 0.1:
        logger.info(f"No significant difference in score distribution between groups with chip_evidence=True and False, p-value={p_value:.3f}")
    else:
        logger.warning("There is a significant difference in score distribution between groups.")
    return balanced_df




def pearsonr_tolerating_nan(x,y):
    x,y=remove_nan_inf(x,y)
    corr,pval=pearsonr(x,y)
    return corr,pval
    
    

def replace_inf(x):
    """
    Replace inf with second maximum unique value of the column * 1.5
    """
    return x.replace(np.inf,x.sort_values(ascending=False).unique()[1]*1.5)
    



def get_minimum_positive(x):
    """
    return the minimum positive value of column x
    """
    return x[x>0].sort_values().unique()[0]
    
    













#------------------------------------------
# Functions for gnomAD enrichment analysis
#------------------------------------------
def bin_and_label(df, 
                  column_name, 
                  bin_edges, 
                  new_column_name="Bin"):
    bin_edges = sorted(set(bin_edges))
    labels = [f"{bin_edges[i]} - {bin_edges[i+1]}" for i in range(len(bin_edges)-1)]
    df[new_column_name] = pd.cut(df[column_name], bins=bin_edges, labels=labels, include_lowest=True)
    return df


def bin_two_columns(df,x_bin_colname,x_bins,
                      hue_bin_colname,hue_bin_item_name_dict,hue_bins,hue_alias):
    """
    Bin two columns in order to calculate odds ratio
    """
    df=bin_and_label(df, x_bin_colname, x_bins, "X_Bin")
    df=bin_and_label(df, hue_bin_colname, hue_bins, hue_alias)
    df[hue_alias]=df[hue_alias].apply(lambda x: hue_bin_item_name_dict[x])
    df=df.loc[:,["X_Bin",hue_alias]].copy().groupby(["X_Bin",hue_alias]).size().unstack()
    df.fillna(0,inplace=True)
    return df


# Odds ratio calculation: in/out bin v.s is/isn't in column
#           common    non-common
# in-bin     A           B
# out-bin    C           D
# 
# odds_common = AD/BC
# A, B may be 0
# C, D can never be 0 
def calc_odds_ratio(df,row_name,col_name):
    A=df.loc[row_name,col_name]
    B=df.loc[row_name,:].sum()-A
    C=df.loc[:,col_name].sum()-A
    D=df.loc[:,].sum().sum()-A-B-C
    table = np.array([[A, B], [C, D]])
    odds_ratio, p = fisher_exact(table)
    se_log_or=np.sqrt(1/A+1/B+1/C+1/D)
    ci_lower=np.exp(np.log(odds_ratio)-1.96*se_log_or)
    ci_upper=np.exp(np.log(odds_ratio)+1.96*se_log_or)
    return odds_ratio, p, ci_lower, ci_upper



def transform_data(df,x_var,hue_var):
    df[x_var]=df.index
    df=df.reset_index(drop=True)
    # melt columns
    or_columns=[i for i in df.columns if i.startswith("or")]
    pval_columns=[i for i in df.columns if i.startswith("pval")]
    ci_low_columns=[i for i in df.columns if i.startswith("ci_low")]
    ci_high_columns=[i for i in df.columns if i.startswith("ci_high")]
    df_or=pd.melt(df,id_vars=[x_var],value_vars=or_columns,var_name=hue_var,value_name="odds_ratio")
    df_pval=pd.melt(df,id_vars=[x_var],value_vars=pval_columns,var_name=hue_var,value_name="pval")
    df_ci_low=pd.melt(df,id_vars=[x_var],value_vars=ci_low_columns,var_name=hue_var,value_name="ci_low")
    df_ci_high=pd.melt(df,id_vars=[x_var],value_vars=ci_high_columns,var_name=hue_var,value_name="ci_high")
    # rename columns
    df_or[hue_var]=df_or[hue_var].str.replace("or_","")
    df_pval[hue_var]=df_pval[hue_var].str.replace("pval_","")
    df_ci_low[hue_var]=df_ci_low[hue_var].str.replace("ci_low_","")
    df_ci_high[hue_var]=df_ci_high[hue_var].str.replace("ci_high_","")
    # merge data frames
    df=pd.merge(df_or,df_pval,on=[x_var,hue_var])
    df=pd.merge(df,df_ci_low,on=[x_var,hue_var])
    df=pd.merge(df,df_ci_high,on=[x_var,hue_var])
    return df



def calc_or(df_orig,x_var,hue_var,out_group=None):
    """
    df: each row is a bin for x axis, each column is a category
    For each grid of df, calculate odds ratio, p-value, confidence interval lower and upper bounds
    """
    # remove rows with any zero
    df=df_orig.loc[(df_orig!=0).all(axis=1),:].copy()    
    categories=df.columns.tolist()
    if out_group is not None and out_group in categories:
        categories.remove(out_group)
    for category in categories:
        or_list = []
        pval_list = []
        ci_low_list = []
        ci_high_list = []
        for index, _ in df.iterrows():
            or_value, pval_value, ci_low_value, ci_high_value = calc_odds_ratio(df_orig, index, category)
            or_list.append(or_value)
            pval_list.append(pval_value)
            ci_low_list.append(ci_low_value)
            ci_high_list.append(ci_high_value)
        df.loc[:,f'or_{category}'] = or_list
        df.loc[:,f'pval_{category}'] = pval_list
        df.loc[:,f'ci_low_{category}'] = ci_low_list
        df.loc[:,f'ci_high_{category}'] = ci_high_list
    df=transform_data(df,x_var,hue_var)
    return df



def plot_or(df_plot, x_colname, y_colname, hue_colname, title, color_mapping, out_name=None, ax=None):
    # determine transparency
    df_plot["alphas"] = df_plot["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))
    if ax is None:
        plt.figure(figsize=(5, 5))
        ax = plt.gca()
    for variant, df_subset in df_plot.groupby(hue_colname):
        if variant not in color_mapping.keys():
            continue
        ax.plot(df_subset[x_colname], df_subset[y_colname], '--', color=color_mapping[variant], label=variant)
        ax.scatter(df_subset[x_colname], df_subset[y_colname], color=color_mapping[variant], alpha=df_subset['alphas'])
        ax.errorbar(df_subset[x_colname], 
                    df_subset[y_colname], 
                    yerr=[df_subset[y_colname]-df_subset['ci_low'],df_subset['ci_high']-df_subset['odds_ratio']],
                    color=color_mapping[variant],
                    capsize=3, markeredgewidth=1)
    ax.set_title(title)
    ax.axhline(y=1, color='black', linestyle=':')
    ax.legend()
    ax.label_outer()  # Only show outer labels to share axis labels
    if out_name is not None:
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(out_name)
        plt.close()




def _convert_interval_to_midpoint(interval_str):
    start, end = map(float, interval_str.split(' - '))
    midpoint = (start + end) / 2
    return midpoint

def plot_or_jitter(df_plot,x_colname,y_colname,title,out_name,color_mapping, jitter_strength=0.08):
    # determine transparency
    df_plot["alphas"]=df_plot["pval"].apply(lambda x: 1.0 if x < 0.001 else (0.1 if x > 0.05 else 0.5))
    plt.figure(figsize=(5, 5)) 
    original_x_labels = {}
    for variant, df_subset in df_plot.groupby('variant_type'):
        if variant not in color_mapping.keys():
            continue
        x_values = df_subset[x_colname]
        y_values = df_subset[y_colname]
        x_values_numeric = np.float_(x_values.apply(_convert_interval_to_midpoint))
        # replace infinity
        diff=x_values_numeric[2]-x_values_numeric[1]
        if np.isinf(x_values_numeric[0]):
            x_values_numeric[0] = x_values_numeric[1]-diff
        if np.isinf(x_values_numeric[-1]):
            x_values_numeric[-1] = x_values_numeric[-2]+diff
        logger.info(x_values_numeric)
        for original_label, midpoint in zip(x_values, x_values_numeric):
            original_x_labels[midpoint] = original_label
        jittered_x_values = x_values_numeric + jitter_strength 
        jitter_strength+=0.08
        
        plt.plot(jittered_x_values, y_values, '--', color=color_mapping[variant], label=variant)
        plt.scatter(jittered_x_values, y_values, color=color_mapping[variant], alpha=df_subset['alphas'])
        plt.errorbar(jittered_x_values, 
                    y_values, 
                    yerr=[y_values - df_subset['ci_low'], df_subset['ci_high'] - y_values],
                    fmt='none', 
                    color=color_mapping[variant],
                    capsize=3, markeredgewidth=1)
    plt.title(title)
    plt.axhline(y=1, color='black', linestyle=':')
    plt.xlabel(x_colname)
    plt.ylabel(y_colname)
    midpoints = sorted(original_x_labels.keys())
    original_labels = [original_x_labels[midpoint] for midpoint in midpoints]
    plt.xticks(midpoints, original_labels, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_name)
    plt.close()