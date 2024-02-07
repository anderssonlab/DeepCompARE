from scipy.spatial import cKDTree
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu
from utils import remove_nan_inf
from scipy.stats import pearsonr


def match_by_score(df, score_col, label_col):
    """
    For each positive sample, find a unique negative sample with the most similar score using a vectorized approach.
    
    Parameters:
    - df: DataFrame containing the data.
    - score_col: Name of the column containing the scores.
    - label_col: Name of the column containing binary labels (1 for positive, 0 for negative).
    
    Returns:
    - DataFrame containing uniquely matched positive and negative samples.
    """
    # Separate positive and negative samples
    logger.info(f"Data frame has shape {df.shape}")
    positives = df[df[label_col] == 1]
    negatives = df[df[label_col] == 0]
    logger.info(f"{len(positives)} positive samples and {len(negatives)} negative samples")
    # Create a KDTree for efficient nearest neighbor search
    tree = cKDTree(negatives[[score_col]])
    # Compute the distance to the nearest neighbor for each positive sample
    _, indices = tree.query(positives[[score_col]], k=1)
    # Ensure unique matches
    unique_negatives = negatives.iloc[np.unique(indices)]
    matched_positives = positives.iloc[np.unique(indices, return_index=True)[1]]
    # Combine matched positive and negative sampless
    matched_df = pd.concat([matched_positives, unique_negatives])
    logger.info(f"Matched data frame has shape {matched_df.shape}")
    return matched_df





def match_by_decile(df,score_col,label_col):
    logger.info(f"Data frame has shape {df.shape}")
    df=df.copy()    
    df['score_decile'] = pd.qcut(df[score_col], 10, labels=False, duplicates='drop')
    
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
        return None
    
    _, p_value = mannwhitneyu(scores_true, scores_false)
    
    if p_value > 0.05:
        logger.info(f"No significant difference in score distribution between groups with chip_evidence=True and False, p-value={p_value:.3f}")
    else:
        logger.warning("There is a significant difference in score distribution between groups.")
    
    return balanced_df




def pearsonr_tolerating_nan(x,y):
    x,y=remove_nan_inf(x,y)
    corr,pval=pearsonr(x,y)
    return corr,pval
    