import pandas as pd
import numpy as np
from loguru import logger
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact
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


#----------------------
# plot odds ratio
#----------------------



def bin_and_label(df, column_name, bin_edges):
    bin_edges = sorted(set(bin_edges))
    labels = [f"{bin_edges[i]} - {bin_edges[i+1]}" for i in range(len(bin_edges)-1)]
    df[f"{column_name}_bin"] = pd.cut(df[column_name], bins=bin_edges, labels=labels, include_lowest=True)
    df[f"{column_name}_bin"] = df[f"{column_name}_bin"].astype(pd.CategoricalDtype(categories=labels, ordered=True))
    return df



def fisher_exact_with_ci(table):
    odds_ratio, p = fisher_exact(table)
    se_log_or=np.sqrt(1/table[0,0]+1/table[0,1]+1/table[1,0]+1/table[1,1])
    ci_lower=np.exp(np.log(odds_ratio)-1.96*se_log_or)
    ci_upper=np.exp(np.log(odds_ratio)+1.96*se_log_or)
    return odds_ratio, p, ci_lower, ci_upper


# Odds ratio calculation: in/out bin v.s True/False
#            True        False
# in-bin     A           B
# out-bin    C           D
# 
# odds_ratio_common = AD/BC
def odds_ratio_one_bin(df,bin_name):
    """
    Calculate odds ratio for a given bin
    """
    A=df.loc[bin_name,"True"]
    B=df.loc[bin_name,"False"]
    C=df.loc[:,"True"].sum()-A
    D=df.loc[:,"False"].sum()-B
    table = np.array([[A, B], [C, D]])
    return fisher_exact_with_ci(table)


def odds_ratio_one_df(df,suffix):
    """
    Iterate over all bins (rows) and calculate odds ratio for each bin
    """
    if df.shape[0]==0:
        logger.warning("Empty data frame")
        return None,None
    df=df.loc[(df["True"]>0) & (df["False"]>0)]
    or_list = []
    pval_list = []
    ci_low_list = []
    ci_high_list = []
    for index, _ in df.iterrows():
        or_value, pval_value, ci_low_value, ci_high_value = odds_ratio_one_bin(df, index)
        or_list.append(or_value)
        pval_list.append(pval_value)
        ci_low_list.append(ci_low_value)
        ci_high_list.append(ci_high_value)
    df_res=pd.DataFrame({f"or":or_list,
                         f"pval":pval_list,
                         f"ci_low":ci_low_list,
                         f"ci_high":ci_high_list},
                        index=df.index)
    df_res["threshold"]=suffix
    return df_res


def calc_or_by_various_thresholds(df,threshold_col,threshold_list,operation,bin_col):
    df_res=pd.DataFrame()
    n_list=[]
    for threshold in threshold_list:
        if operation=="larger":
            df["target"]=(df[threshold_col]>threshold)
        if operation=="smaller":
            df["target"]=(df[threshold_col]<threshold)
        n_list.append(df["target"].sum())
        df["target"]=df["target"].astype(str)
        df_plot=df.groupby([bin_col,"target"]).size().unstack()
        df_temp=odds_ratio_one_df(df_plot,threshold)
        df_res=pd.concat([df_res,df_temp])
    return df_res,n_list




