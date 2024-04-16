import pandas as pd
import numpy as np
from loguru import logger
from scipy.stats import mannwhitneyu
from utils import remove_nan_inf
from scipy.stats import pearsonr, fisher_exact




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
    
    
    
    
def bin_and_label(df, column_name, bin_edges):
    bin_edges = sorted(set(bin_edges))
    labels = [f"{bin_edges[i]}-{bin_edges[i+1]}" for i in range(len(bin_edges)-1)]
    df['Bin'] = pd.cut(df[column_name], bins=bin_edges, labels=labels, include_lowest=True)
    return df






# Odds ratio calculation: in/out bin v.s is/isn't in column
#           common    non-common
# in-bin     A           B
# out-bin    C           D
# odds_rare = AD/BC
# A, B may be 0
# C, D can never be 0 
def calc_odds_ratio(df,row_name,col_name):
    A=df.loc[row_name,col_name]
    B=df.loc[row_name,:].sum()-A
    C=df.loc[:,col_name].sum()-A
    D=df.loc[:,].sum().sum()-A-B-C
    table = np.array([[A, B], [C, D]])
    odds_ratio, p = fisher_exact(table)
    return odds_ratio, p