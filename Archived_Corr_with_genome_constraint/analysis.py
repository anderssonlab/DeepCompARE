"""
For each region in /isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_K562.bed
Get avg z score from /isdata/alab/people/pcr980/Resource/constraint_z_genome_1kb.qc.download.txt
Then get max gradxinp from /isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/gradxinp_CAGE_K562.csv
Then getmax ism from /isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/ism_CAGE_K562.csv
Calculate corr(z, gradxinp) and corr(z, ism)
"""


import pandas as pd
import numpy as np
from loguru import logger
from scipy.stats import pearsonr

import sys
sys.path.insert(1,"/isdata/alab/people/pcr980/DeepCompare/Scripts_python/")
from utils import read_featimp



def get_avg_zscore(region, zscore_df):
    """
    Get avg z score from zscore_df for region
    """
    chrom=region[0]
    start=int(region[1])
    end=int(region[2])
    zscore_df = zscore_df[(zscore_df['chrom'] == chrom) & (zscore_df['end'] >= start) & (zscore_df['start'] <= end)]
    # if zscore is empty, return np.nan
    if zscore_df.empty:
        return np.nan
    zscore_df['overlap_len'] = zscore_df.apply(lambda row: min(end, row['end']) - max(start, row['start']), axis=1)
    zscore_df['weighted_zscore'] = zscore_df['overlap_len'] * zscore_df['z']
    return zscore_df['weighted_zscore'].sum() / (end - start + 1)


if __name__=='__main__':
    logger.info("Start analysis")
    zscore_df = pd.read_csv('/isdata/alab/people/pcr980/Resource/constraint_z_genome_1kb.qc.download.txt', sep='\t')
    logger.info("Read in z score file")
    
    gradxinp_df = read_featimp('/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/gradxinp_CAGE_K562.csv', track_num=1)
    logger.info("Read in gradxinp file")
    ism_df = read_featimp('/isdata/alab/people/pcr980/DeepCompare/Pd3_metadata_and_featimp/ism_CAGE_K562.csv', track_num=1)
    logger.info("Read in ism file")
    
    # decide whether to use abs
    # gradxinp_df = gradxinp_df.abs()
    # ism_df = ism_df.abs()
    
    regions_df = pd.read_csv('/isdata/alab/people/pcr980/DeepCompare/Pd1_bed_processed/CAGE_K562.bed', sep='\t', header=None)
    logger.info("Read in regions file")
    regions_df.columns = ['chrom', 'start', 'end',"x","score","strand"]
    regions_df = regions_df.drop(columns=['x', 'score', 'strand'])
    regions_df['avg_zscore'] = regions_df.apply(lambda row: get_avg_zscore((row[0],row[1],row[2]), zscore_df), axis=1)
    logger.info("Get avg z score for each region")
    
    regions_df['max_gradxinp'] = gradxinp_df.values.max(axis=1)
    logger.info("Get max gradxinp for each region")
    
    regions_df['max_ism'] = ism_df.values.max(axis=1)
    logger.info("Get max ism for each region")
    
    # remove rows containing nan
    regions_df.to_csv('/isdata/alab/people/pcr980/DeepCompare/Corr_with_genome_constraint/regions_df.csv', index=False)
    logger.info(f"Before remove nan: {regions_df.shape}")
    regions_df = regions_df.dropna()
    logger.info(f"After remove nan: {regions_df.shape}")
    
    corr, pval = pearsonr(regions_df['avg_zscore'], regions_df['max_gradxinp'])
    logger.info(f"Correlation between z score and gradxinp: {corr}, p-value: {pval}")
    corr, pval = pearsonr(regions_df['avg_zscore'], regions_df['max_ism'])
    logger.info(f"Correlation between z score and ism: {corr}, p-value: {pval}")
    