import pandas as pd
import numpy as np
import pyBigWig

from loguru import logger


# TODO: add true CAGE/DHS/STARR/CAGE as well

def calculate_average_bw(bw, region):
    chrom, start, end = region
    values = bw.values(chrom, start, end, numpy=True)
    # calculate averge signal on non-NaN values
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0


def annotate(df,cell_line,out_path):
    for histone_modification in ['H2AFZ', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9ac', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K79me2', 'H4K20me1']:
        bw = pyBigWig.open(f"Bigwig_data/{histone_modification}_{cell_line}.bigwig")
        df[f'log_signal_{histone_modification}'] = df.apply(lambda row: calculate_average_bw(bw, (row[0], row[1], row[2])), axis=1)
        logger.info(f'Done with {histone_modification}')
    df.to_csv(out_path,index=False)


for re in ["promoters", "enhancers"]:
    for cell_line in ["hepg2","k562"]:
        df=pd.read_csv(f'/isdata/alab/people/pcr980/DeepCompare/Pd4_promoters_enhancers_and_featimp/{re}_{cell_line}.bed', sep='\t', header=None)
        # retain first 3 columns
        df = df.iloc[:, :3]
        df.columns = ['chrom', 'start', 'end']
        annotate(df,cell_line,f'{re}_{cell_line}.csv')
        
        

# nohup python3 quantify_promoters_enhancers.py > quantify_promoters_enhancers.log &