import pandas as pd
import numpy as np
import pyBigWig
import os

from loguru import logger


def calculate_average_bw(bw, region):
    chrom, start, end = region
    values = bw.values(chrom, start, end, numpy=True)
    # calculate averge signal on non-NaN values
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0


def annotate(df,out_path):
    chip_path="/isdata/alab/people/pcr980/Resource/ChIP_Seq/"
    # get all file names ending with .bigWig
    chip_list = [f.split('_')[0] for f in os.listdir(chip_path) if f.endswith('.bigWig')]
    # remove the .bigWig extension
    chip_list = [f.split('.')[0] for f in chip_list]
    for tf_chip in chip_list:
        bw = pyBigWig.open(f"{chip_path}{tf_chip}.bigWig")
        df[f'chip_{tf_chip}'] = df.apply(lambda row: calculate_average_bw(bw, (row[0], row[1], row[2])), axis=1)
        logger.info(f'Done with {tf_chip}')
    df.to_csv(out_path,index=False)




for re in ["proximal", "distal"]:
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/dhs_{re}_k562.tsv', sep=' ')
    # retain first 3 columns
    df = df.iloc[:, :3]
    df.columns = ['chrom', 'start', 'end']
    annotate(df,f'chip_{re}.csv')
    logger.info(f'Done with {re}')
    



for re in ["promoters", "enhancers"]:
    df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{re}_k562.bed', sep='\t', header=None)
    # retain first 3 columns
    df = df.iloc[:, :3]
    df.columns = ['chrom', 'start', 'end']
    annotate(df,f'chip_{re}.csv')
    logger.info(f'Done with {re}')
    
    

# nohup python3 quantify_chip.py > quantify_chip.log &